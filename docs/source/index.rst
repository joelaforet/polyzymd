PolyzyMD
========

.. container:: pz-hero

   .. image:: _static/polyzymd_badge_256.png
      :alt: PolyzyMD logo: an enzyme surrounded by polymer chains
      :class: dark-light

   .. container::

      **Build, run and analyze molecular dynamics of proteins** — alone or with
      polymers, ligands and co-solvents — on OpenMM or GROMACS, with every
      build checked against its config and every analysis compared across
      replicates.

      .. container:: pz-buttons

         .. button-ref:: get_started/installation
            :color: primary
            :shadow:

            Install PolyzyMD

         .. button-ref:: get_started/quickstart
            :color: primary
            :outline:

            Run your first simulation

How PolyzyMD works
------------------

.. raw:: html

   <nav class="pz-flow" aria-label="The PolyzyMD workflow">
     <a class="pz-step" href="tutorials/prepare_pdb_for_openff.html"><span class="pz-step-n">1</span><span class="pz-step-t">Prepare</span><span class="pz-step-d">a clean PDB, monomer SMILES</span></a>
     <a class="pz-step" href="reference/configuration.html"><span class="pz-step-n">2</span><span class="pz-step-t">Configure</span><span class="pz-step-d">one config.yaml per condition</span></a>
     <a class="pz-step" href="explanation/physics_audit.html"><span class="pz-step-n">3</span><span class="pz-step-t">Build &amp; audit</span><span class="pz-step-d">every term checked against the config</span></a>
     <a class="pz-step" href="how_to/hpc_slurm.html"><span class="pz-step-n">4</span><span class="pz-step-t">Run</span><span class="pz-step-d">OpenMM or GROMACS, laptop or SLURM</span></a>
     <a class="pz-step" href="how_to/analysis_chooser.html"><span class="pz-step-n">5</span><span class="pz-step-t">Analyze</span><span class="pz-step-d">per replicate, with uncertainty</span></a>
     <a class="pz-step" href="how_to/study_freeze.html"><span class="pz-step-n">6</span><span class="pz-step-t">Compare &amp; publish</span><span class="pz-step-d">conditions, statistics, a frozen deposit</span></a>
   </nav>

Find what you need
------------------

The documentation is organized by what you want to do right now.

.. grid:: 1 1 2 2
   :gutter: 3

   .. grid-item-card:: :octicon:`mortar-board;1.2em` Tutorials
      :link: tutorials/index
      :link-type: doc

      Learn by doing: guided lessons from a PDB file to a compared study.

      +++
      For when you are new to PolyzyMD

   .. grid-item-card:: :octicon:`tools;1.2em` How-to guides
      :link: how_to/index
      :link-type: doc

      Recipes for one task: polymers, restraints, SLURM, GROMACS, running
      longer, analyses, publishing a study.

      +++
      For when you know what you want to do

   .. grid-item-card:: :octicon:`light-bulb;1.2em` Explanation
      :link: explanation/index
      :link-type: doc

      Why PolyzyMD works the way it does: the physics audit, replicates and
      statistics, study design.

      +++
      For when you want to understand

   .. grid-item-card:: :octicon:`book;1.2em` Reference
      :link: reference/index
      :link-type: doc

      Every command, config key, file, check and analysis function, and the
      papers behind them.

      +++
      For when you need a fact

Common tasks
------------

.. container:: pz-buttons

   .. button-ref:: get_started/examples
      :color: primary
      :outline:

      Run the examples

   .. button-ref:: how_to/dynamic_polymers
      :color: primary
      :outline:

      Build polymers from SMILES

   .. button-ref:: how_to/cosolvents
      :color: primary
      :outline:

      Add a co-solvent

   .. button-ref:: how_to/run_gromacs
      :color: primary
      :outline:

      Run with GROMACS

   .. button-ref:: local-nvidia-gpu
      :ref-type: ref
      :color: primary
      :outline:

      Run on your NVIDIA GPU

   .. button-ref:: how_to/hpc_slurm
      :color: primary
      :outline:

      Run on a SLURM cluster

   .. button-ref:: how_to/analysis_compare_conditions
      :color: primary
      :outline:

      Compare conditions

   .. button-ref:: how_to/own_analysis
      :color: primary
      :outline:

      Use your own MDAnalysis scripts

   .. button-ref:: contributor_guide/index
      :color: primary
      :outline:

      Contribute to PolyzyMD

New in 1.3
----------

.. grid:: 1 1 3 3
   :gutter: 2

   .. grid-item-card:: Physics audit
      :link: explanation/physics_audit
      :link-type: doc

      Every build is compared, term by term, with the force fields and charges
      its config names.

   .. grid-item-card:: Polymers from SMILES
      :link: tutorials/star_sited_monomers
      :link-type: doc

      Star-sited monomers, drawn end caps and tacticity, built with mBuild.

   .. grid-item-card:: Stop, continue, run longer
      :link: how_to/run_longer
      :link-type: doc

      Interrupted runs continue on OpenMM and GROMACS; finished OpenMM runs
      extend to a longer production.

.. note::

   PolyzyMD is under active development. Found a problem?
   `Open an issue on GitHub <https://github.com/joelaforet/polyzymd/issues>`_.

.. toctree::
   :hidden:
   :maxdepth: 2

   Get started <get_started/index>
   Tutorials <tutorials/index>
   How-to guides <how_to/index>
   Reference <reference/index>
   Explanation <explanation/index>
   Contributing <contributor_guide/index>
