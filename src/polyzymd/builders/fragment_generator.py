"""
Fragment generator for dynamic polymer generation.

This module handles the generation of polymer fragments from raw monomer SMILES
using ATRP (Atom-Transfer Radical Polymerization) reaction chemistry.

The workflow mirrors the experimental process:
1. Load raw monomer SMILES (unactivated)
2. Run initiation reaction (chlorination for ATRP)
3. Run polymerization reaction to generate fragments
4. Terminate 1-site fragments (restore alkene)
5. Build and cache MonomerGroup for polymer building

Made by PolyzyMD, by Joseph R. Laforet Jr.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import rdFMCS

from polymerist.rdutils.reactions.reactions import AnnotatedReaction
from polymerist.rdutils.reactions.reactors import PolymerizationReactor
from polymerist.rdutils.sanitization import explicit_mol_from_SMILES, Smiles
from polymerist.polymers.monomers import MonomerGroup
from polymerist.rdutils.bonding.portlib import get_num_ports
from polymerist.genutils.iteration import sort_dict_by_keys

logger = logging.getLogger(__name__)


class FragmentGenerator:
    """Generates polymer fragments from raw monomer SMILES using ATRP chemistry.

    This class encapsulates the fragment generation workflow, which involves:
    1. Loading reaction templates from .rxn files
    2. Running initiation reactions on raw monomers
    3. Generating all possible fragments via polymerization reactions
    4. Terminating 1-site fragments to restore terminal alkenes
    5. Creating a MonomerGroup ready for polymer building

    Attributes:
        initiation_rxn_path: Path to ATRP initiation reaction file
        polymerization_rxn_path: Path to ATRP polymerization reaction file
        termination_rxn_path: Path to ATRP termination reaction file
        cache_directory: Directory for caching generated fragments
    """

    def __init__(
        self,
        initiation_rxn_path: Path,
        polymerization_rxn_path: Path,
        termination_rxn_path: Path,
        cache_directory: Path,
    ):
        """Initialize the fragment generator.

        Args:
            initiation_rxn_path: Path to ATRP initiation .rxn file
            polymerization_rxn_path: Path to ATRP polymerization .rxn file
            termination_rxn_path: Path to ATRP termination .rxn file
            cache_directory: Directory for caching MonomerGroup JSON
        """
        self.initiation_rxn_path = Path(initiation_rxn_path)
        self.polymerization_rxn_path = Path(polymerization_rxn_path)
        self.termination_rxn_path = Path(termination_rxn_path)
        self.cache_directory = Path(cache_directory)

        # Lazy-loaded reactions
        self._initiation_rxn: Optional[AnnotatedReaction] = None
        self._polymerization_rxn: Optional[AnnotatedReaction] = None
        self._termination_rxn: Optional[AnnotatedReaction] = None

        # Validate paths exist
        self._validate_paths()

    def _validate_paths(self) -> None:
        """Validate that all reaction files exist."""
        for path, name in [
            (self.initiation_rxn_path, "initiation"),
            (self.polymerization_rxn_path, "polymerization"),
            (self.termination_rxn_path, "termination"),
        ]:
            if not path.exists():
                raise FileNotFoundError(f"ATRP {name} reaction file not found: {path}")

    @property
    def initiation_rxn(self) -> AnnotatedReaction:
        """Load initiation reaction lazily."""
        if self._initiation_rxn is None:
            self._initiation_rxn = AnnotatedReaction.from_rxnfile(self.initiation_rxn_path)
            logger.debug(f"Loaded initiation reaction from {self.initiation_rxn_path}")
        return self._initiation_rxn

    @property
    def polymerization_rxn(self) -> AnnotatedReaction:
        """Load polymerization reaction lazily."""
        if self._polymerization_rxn is None:
            self._polymerization_rxn = AnnotatedReaction.from_rxnfile(self.polymerization_rxn_path)
            logger.debug(f"Loaded polymerization reaction from {self.polymerization_rxn_path}")
        return self._polymerization_rxn

    @property
    def termination_rxn(self) -> AnnotatedReaction:
        """Load termination reaction lazily."""
        if self._termination_rxn is None:
            self._termination_rxn = AnnotatedReaction.from_rxnfile(self.termination_rxn_path)
            logger.debug(f"Loaded termination reaction from {self.termination_rxn_path}")
        return self._termination_rxn

    def _partial_substruct_match(
        self, mol: Chem.Mol, query: Chem.Mol, threshold: float = 0.8
    ) -> bool:
        """Check if mol partially matches query using MCS (Maximum Common Substructure).

        Args:
            mol: Molecule to check
            query: Query molecule (fragment)
            threshold: Minimum fraction of query atoms that must match

        Returns:
            True if match fraction >= threshold
        """
        mcs = rdFMCS.FindMCS(
            [mol, query],
            completeRingsOnly=False,
            matchValences=False,
            ringMatchesRingOnly=False,
            timeout=5,
        )
        if not mcs.smartsString:
            return False
        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
        n_mcs_atoms = mcs_mol.GetNumAtoms()
        n_query_atoms = query.GetNumAtoms()
        return (n_mcs_atoms / n_query_atoms) >= threshold

    def _find_parent_monomer(
        self,
        monomers: Dict[str, Chem.Mol],
        query_mol: Chem.Mol,
        threshold: float = 0.77,
    ) -> Optional[str]:
        """Find the parent monomer that contains a given fragment.

        Args:
            monomers: Dictionary of monomer name -> RDKit Mol
            query_mol: Fragment molecule to match
            threshold: Minimum match threshold

        Returns:
            Name of the parent monomer, or None if not found
        """
        for name, monomer in monomers.items():
            if self._partial_substruct_match(monomer, query_mol, threshold=threshold):
                return name
        return None

    def generate_fragments(
        self,
        monomer_smiles: Dict[str, str],
        polymerization_depth: int = 5,
        match_threshold: float = 0.77,
    ) -> MonomerGroup:
        """Generate polymer fragments from raw monomer SMILES.

        This method executes the full ATRP fragment generation workflow:
        1. Parse raw monomer SMILES into RDKit molecules
        2. Add Cl2 for initiation reaction
        3. Run initiation to activate monomers
        4. Run polymerization to generate all fragment types
        5. Name fragments by parent monomer and functionality (1-site, 2-site)
        6. Terminate 1-site fragments to restore alkene
        7. Create final MonomerGroup

        Args:
            monomer_smiles: Dictionary of monomer name -> SMILES string
            polymerization_depth: Maximum reaction depth for fragment enumeration
            match_threshold: Threshold for matching fragments to parent monomers

        Returns:
            MonomerGroup containing all named fragments ready for polymer building
        """
        logger.info(f"Generating fragments for {len(monomer_smiles)} monomers")

        # Step 1: Convert SMILES to RDKit molecules
        monomers: Dict[str, Chem.Mol] = {
            name: explicit_mol_from_SMILES(smiles) for name, smiles in monomer_smiles.items()
        }

        # Add Cl2 for initiation
        monomers["Cl2"] = explicit_mol_from_SMILES("[Cl]-[Cl]")

        # Step 2: Run initiation reaction
        logger.info("Running initiation reactions...")
        initiation_reactor = PolymerizationReactor(self.initiation_rxn)

        # Get starting monomers (before initiation)
        starting_monomers = initiation_reactor.propagate_pooled(
            monomers.values(),
            rxn_depth_max=0,
            clear_dummy_labels=True,
        )

        # Get activated monomers (after initiation)
        activated_monomers = initiation_reactor.propagate_pooled(
            monomers.values(),
            rxn_depth_max=1,
            clear_dummy_labels=True,
        )

        # Remove starting monomers from activated set
        for smiles in list(starting_monomers.keys()):
            activated_monomers.pop(smiles, None)

        logger.info(f"Generated {len(activated_monomers)} activated monomers")

        # Step 3: Name activated monomers by parent
        activated_named_monomers: Dict[str, Chem.Mol] = {}
        for canon_smiles, fragment_mol in activated_monomers.items():
            query_mol = Chem.MolFromSmarts(canon_smiles)
            parent_name = self._find_parent_monomer(monomers, query_mol, threshold=match_threshold)
            if parent_name and parent_name != "Cl2":
                activated_named_monomers[parent_name] = fragment_mol
                logger.debug(f"Matched activated monomer to parent: {parent_name}")

        # Step 4: Run polymerization reactions
        logger.info("Running polymerization reactions...")
        polymerization_reactor = PolymerizationReactor(self.polymerization_rxn)
        activated_fragments = polymerization_reactor.propagate_pooled(
            activated_named_monomers.values(),
            rxn_depth_max=polymerization_depth,
            clear_dummy_labels=True,
            allow_resampling=True,
        )

        logger.info(f"Generated {len(activated_fragments)} fragment types")

        # Step 5: Name fragments by parent and functionality
        from string import ascii_lowercase
        from collections import defaultdict, Counter

        named_fragments: Dict[str, Smiles] = {}
        fragment_name_modifiers: Dict[str, Counter] = defaultdict(Counter)

        for canon_smiles, fragment_mol in activated_fragments.items():
            functionality = get_num_ports(fragment_mol)
            if functionality == 0:
                continue  # Skip fragments with no ports

            query_mol = Chem.MolFromSmarts(canon_smiles)
            parent_name = self._find_parent_monomer(monomers, query_mol, threshold=match_threshold)
            if parent_name is None or parent_name == "Cl2":
                continue

            # Generate unique suffix for this (parent, functionality) pair
            suffix = ascii_lowercase[fragment_name_modifiers[parent_name][functionality]]
            fragment_name = f"{parent_name}-{functionality}{suffix}"
            named_fragments[fragment_name] = canon_smiles
            fragment_name_modifiers[parent_name][functionality] += 1
            logger.debug(f"Named fragment: {fragment_name}")

        # Step 6: Terminate 1-site fragments
        logger.info("Running termination reactions on 1-site fragments...")
        terminated_fragments = self._terminate_one_site_fragments(
            named_fragments, activated_fragments, monomers, match_threshold
        )

        # Step 7: Build final fragment dictionary
        final_fragments = self._build_final_fragments(
            named_fragments, terminated_fragments, monomers, match_threshold
        )

        # Create MonomerGroup
        monogrp = MonomerGroup(sort_dict_by_keys(final_fragments, reverse=True))
        logger.info(f"Created MonomerGroup with {len(final_fragments)} named fragments")

        return monogrp

    def _terminate_one_site_fragments(
        self,
        named_fragments: Dict[str, Smiles],
        activated_fragments: Dict[str, Chem.Mol],
        monomers: Dict[str, Chem.Mol],
        match_threshold: float,
    ) -> Dict[str, Chem.Mol]:
        """Terminate 1-site fragments to restore the terminal alkene.

        Args:
            named_fragments: Dictionary of fragment name -> SMILES
            activated_fragments: Dictionary of SMILES -> RDKit Mol
            monomers: Original monomer dictionary
            match_threshold: Threshold for parent matching

        Returns:
            Dictionary of terminated fragment SMILES -> RDKit Mol
        """
        # Build activated named fragments dict for termination
        activated_named_fragments: Dict[str, Chem.Mol] = {
            name: explicit_mol_from_SMILES(smiles) for name, smiles in named_fragments.items()
        }

        # Add HCl for termination reaction
        activated_named_fragments["HCl"] = Chem.MolFromSmarts("[H]-[Cl]")

        # Run termination
        termination_reactor = PolymerizationReactor(self.termination_rxn)
        terminated_fragments = termination_reactor.propagate_pooled(
            activated_named_fragments.values(),
            clear_dummy_labels=True,
        )

        # Remove non-terminated fragments
        for smiles in list(activated_fragments.keys()):
            terminated_fragments.pop(smiles, None)
        terminated_fragments.pop("Cl", None)

        return terminated_fragments

    def _build_final_fragments(
        self,
        named_fragments: Dict[str, Smiles],
        terminated_fragments: Dict[str, Chem.Mol],
        monomers: Dict[str, Chem.Mol],
        match_threshold: float,
    ) -> Dict[str, Smiles]:
        """Build the final fragment dictionary with terminated 1-sites.

        For 1-site fragments, use the terminated versions.
        For 2-site fragments, use the original activated versions.

        Args:
            named_fragments: Original named fragments
            terminated_fragments: Terminated fragment molecules
            monomers: Original monomer dictionary
            match_threshold: Threshold for parent matching

        Returns:
            Final dictionary of fragment name -> SMILES
        """
        final_fragments: Dict[str, Smiles] = {}

        # Add 2-site fragments as-is
        for name, smiles in named_fragments.items():
            if "-2" in name:  # 2-site fragment
                final_fragments[name] = smiles

        # Add terminated 1-site fragments
        for canon_smiles, fragment_mol in terminated_fragments.items():
            functionality = get_num_ports(fragment_mol)
            if functionality != 1:
                continue

            query_mol = Chem.MolFromSmarts(canon_smiles)
            parent_name = self._find_parent_monomer(monomers, query_mol, threshold=match_threshold)
            if parent_name and parent_name != "Cl2":
                fragment_name = f"{parent_name}-1a"  # Terminated 1-site
                final_fragments[fragment_name] = canon_smiles

        return final_fragments

    def get_cache_path(self, type_prefix: str) -> Path:
        """Get the path for caching MonomerGroup JSON.

        Args:
            type_prefix: Polymer type prefix for naming

        Returns:
            Path to the MonomerGroup JSON cache file
        """
        self.cache_directory.mkdir(parents=True, exist_ok=True)
        return self.cache_directory / f"{type_prefix}_monomer_group.json"

    def load_or_generate(
        self,
        monomer_smiles: Dict[str, str],
        type_prefix: str,
        force_regenerate: bool = False,
    ) -> MonomerGroup:
        """Load cached MonomerGroup or generate if not cached.

        Args:
            monomer_smiles: Dictionary of monomer name -> SMILES
            type_prefix: Polymer type prefix for cache naming
            force_regenerate: If True, regenerate even if cache exists

        Returns:
            MonomerGroup ready for polymer building
        """
        cache_path = self.get_cache_path(type_prefix)

        if cache_path.exists() and not force_regenerate:
            logger.info(f"Loading cached MonomerGroup from {cache_path}")
            return MonomerGroup.from_file(cache_path)

        # Generate new MonomerGroup
        monogrp = self.generate_fragments(monomer_smiles)

        # Cache to JSON
        monogrp.to_file(cache_path)
        logger.info(f"Cached MonomerGroup to {cache_path}")

        return monogrp
