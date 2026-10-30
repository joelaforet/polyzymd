# Methods and references

This list gives the published source of each method that the PolyzyMD
analyses implement or adapt. An entry marked "not implemented" names a method
that PolyzyMD does not run. Do not cite it as a method of PolyzyMD.

```{note}
For contributors: when you add a method, add its reference here and in the
`References` section of the module docstring.
```

(cite-polyzymd)=
## Cite PolyzyMD

```{index} !cite polyzymd, !how to cite
```

Cite PolyzyMD with the version you used. `polyzymd --version` prints it, and
the "Methods to cite" list of each `build_provenance.md` starts with PolyzyMD
at the version that made the build.

> Laforet Jr., Joseph R., PolyzyMD: build, run and analyze molecular dynamics
> of complex biomolecular systems (version 1.3.0).
> https://github.com/joelaforet/polyzymd. https://doi.org/10.5281/zenodo.23680214

```bibtex
@software{laforet_polyzymd,
  author  = {Laforet Jr., Joseph R.},
  title   = {{PolyzyMD}: build, run and analyze molecular dynamics of complex
             biomolecular systems},
  version = {1.3.0},
  year    = {2026},
  doi     = {10.5281/zenodo.23680214},
  url     = {https://github.com/joelaforet/polyzymd}
}
```

`10.5281/zenodo.23680214` is PolyzyMD's concept DOI on Zenodo: it always resolves to the
newest release. Each release also has a DOI of its own, listed on its Zenodo
page; version 1.3.0 is
[`10.5281/zenodo.23680215`](https://doi.org/10.5281/zenodo.23680215).

The same entry is in
[`CITATION.cff`](https://github.com/joelaforet/polyzymd/blob/main/CITATION.cff),
from which GitHub's "Cite this repository" button makes APA and BibTeX. Cite
also the methods that your build and your analyses used, listed below.

## Building the system

Cite each of these methods that your build used. `build_provenance.md` lists
the ones a build used under "Methods to cite", together with the force fields,
the water model, the ion parameters and the charge models below. For contributors:
each reference and when a build cites it is one `Reference` entry of the
`REFERENCES` table in `polyzymd/simulation/build_citations.py`; to cite one
more method, add an entry there and the reference here.

- Eastman, P., Galvelis, R., Peláez, R. P., Abreu, C. R. A., Farr, S. E.,
  Gallicchio, E., Gorenko, A., Henry, M. M., Hu, F., Huang, J., Krämer, A.,
  Michel, J., Mitchell, J. A., Pande, V. S., Rodrigues, J. P. G. L. M.,
  Rodriguez-Guerra, J., Simmonett, A. C., Singh, S., Swails, J., Turner, P.,
  Wang, Y., Zhang, I., Chodera, J. D., De Fabritiis, G. and Markland, T. E.
  (2024). OpenMM 8: molecular dynamics simulation with machine learning
  potentials. *Journal of Physical Chemistry B*, 128(1), 109-116.
  doi:10.1021/acs.jpcb.3c06662. For a build for OpenMM (the default engine).
- Abraham, M. J., Murtola, T., Schulz, R., Páll, S., Smith, J. C., Hess, B.
  and Lindahl, E. (2015). GROMACS: high performance molecular simulations
  through multi-level parallelism from laptops to supercomputers. *SoftwareX*,
  1-2, 19-25. doi:10.1016/j.softx.2015.06.001. For a build for GROMACS.
- Every build uses the OpenFF Toolkit and OpenFF Interchange to assign
  parameters. They are cited as OpenFF asks
  (<https://openforcefield.org/science/how-to-cite/>): the Toolkit with
  Mobley, D. L., Bannan, C. C., Rizzi, A., Bayly, C. I., Chodera, J. D.,
  Lim, V. T., Lim, N. M., Beauchamp, K. A., Slochower, D. R., Shirts, M. R.,
  Gilson, M. K. and Eastman, P. K. (2018). Escaping atom types in force fields
  using direct chemical perception. *Journal of Chemical Theory and
  Computation*, 14(11), 6076-6092. doi:10.1021/acs.jctc.8b00640, and its
  release, doi:10.5281/zenodo.597754; Interchange, which has no paper, with
  its release, doi:10.5281/zenodo.8147764; and the Open Force Field
  Initiative with Wang, L., Behara, P. K., Thompson, M. W., Gokey, T.,
  Wang, Y., Wagner, J. R., Cole, D. J., Gilson, M. K., Shirts, M. R. and
  Mobley, D. L. (2024). The Open Force Field Initiative: open software and
  open science for molecular modeling. *Journal of Physical Chemistry B*,
  128(29), 7043-7067. doi:10.1021/acs.jpcb.4c01558. An OpenFF Sage or
  Parsley file is cited with its paper and the openff-forcefields release,
  doi:10.5281/zenodo.3473553.

  The three Zenodo DOIs above stand for every release of each package.
  `build_provenance.json` records the versions a build used
  (`openff_toolkit_version`, `openff_interchange_version`,
  `openff_forcefields_version`), and "Methods to cite" gives each with the
  DOI of that release, from a table PolyzyMD ships
  (`polyzymd/data/openff_release_dois.json`, refreshed from Zenodo before
  each PolyzyMD release). A version newer than the table is given with the
  DOI of every release.

- Martínez, L., Andrade, R., Birgin, E. G. and Martínez, J. M. (2009).
  PACKMOL: a package for building initial configurations for molecular
  dynamics simulations. *Journal of Computational Chemistry*, 30(13),
  2157-2164. doi:10.1002/jcc.21224. Used by every build: Packmol places the
  water, co-solvents and polymer chains.
- Machado, M. R. and Pantano, S. (2020). Split the charge difference in two!
  A rule of thumb for adding proper amounts of ions in MD simulations.
  *Journal of Chemical Theory and Computation*, 16(3), 1367-1372.
  doi:10.1021/acs.jctc.9b00953. Used by `solvent.ions.ionic_strength`: the
  neutralizing charge is split between Na+ and Cl-, so the counter-ions count
  toward the ionic strength.
- Klein, C., Sallai, J., Jones, T. J., Iacovella, C. R., McCabe, C. and
  Cummings, P. T. (2016). A hierarchical, component based approach to
  screening properties of soft matter. In *Foundations of Molecular Modeling
  and Simulation* (Molecular Modeling and Simulation), 79-92. Springer
  Singapore. doi:10.1007/978-981-10-1128-3_5. Used by dynamic polymers:
  mBuild builds the coordinates of each chain. PolyzyMD installs mBuild from
  github.com/joelaforet/mbuild (tag `biopolymers-v0.1.0`); the build record
  names its commit.

## Force fields, water, ions and charge models

The build cites each force-field file by its name: ff14SB for
`ff14sb_off_impropers_*.offxml` (the OpenFF port of Amber ff14SB), the Sage
paper for `openff-2.*.offxml`, and the Parsley paper for `openff-1.*.offxml`.
The Sage paper describes 2.0.0; PolyzyMD cites it for every Sage 2 release.
For Sage 2.3.0 (`openff-2.3.0*.offxml`) it also cites the Sage 2.3.0 paper.
For any other file, cite its source.

- Maier, J. A., Martinez, C., Kasavajhala, K., Wickstrom, L., Hauser, K. E.
  and Simmerling, C. (2015). ff14SB: improving the accuracy of protein side
  chain and backbone parameters from ff99SB. *Journal of Chemical Theory and
  Computation*, 11(8), 3696-3713. doi:10.1021/acs.jctc.5b00255
- Boothroyd, S., Behara, P. K., Madin, O. C., Hahn, D. F., Jang, H.,
  Gapsys, V., Wagner, J. R., Horton, J. T., Dotson, D. L., Thompson, M. W.,
  Maat, J., Gokey, T., Wang, L.-P., Cole, D. J., Gilson, M. K., Chodera, J. D.,
  Bayly, C. I., Shirts, M. R. and Mobley, D. L. (2023). Development and
  benchmarking of Open Force Field 2.0.0: the Sage small molecule force field.
  *Journal of Chemical Theory and Computation*, 19(11), 3251-3275.
  doi:10.1021/acs.jctc.3c00039
- Wang, L., Alibay, I., Behara, P. K., Boothroyd, S., Cavender, C. E.,
  Horton, J. T., McIsaac, A. R., Mitchell, A., Morales, B., Thompson, M. W.,
  Wagner, J. R., Westbrook, B. R., Bayly, C. I., Chodera, J. D., Cole, D. J.,
  Eastwood, J. R. B., Shirts, M. R. and Mobley, D. L. (2026). Developing and
  benchmarking Sage 2.3.0 with the AshGC neural network charge model.
  *Journal of Chemical Theory and Computation*, 22(9), 4507-4531.
  doi:10.1021/acs.jctc.6c00169. For `openff-2.3.0*.offxml`, and for the
  AshGC charge model (see below).
- Qiu, Y., Smith, D. G. A., Boothroyd, S., Jang, H., Hahn, D. F., Wagner, J.,
  Bannan, C. C., Gokey, T., Lim, V. T., Stern, C. D., Rizzi, A., Tjanaka, B.,
  Tresadern, G., Lucas, X., Shirts, M. R., Gilson, M. K., Chodera, J. D.,
  Bayly, C. I., Mobley, D. L. and Wang, L.-P. (2021). Development and
  benchmarking of Open Force Field v1.0.0, the Parsley small-molecule force
  field. *Journal of Chemical Theory and Computation*, 17(10), 6262-6280.
  doi:10.1021/acs.jctc.1c00571
- Jorgensen, W. L., Chandrasekhar, J., Madura, J. D., Impey, R. W. and
  Klein, M. L. (1983). Comparison of simple potential functions for simulating
  liquid water. *Journal of Chemical Physics*, 79(2), 926-935.
  doi:10.1063/1.445869. For `solvent.water.model: tip3p`.
- Izadi, S. and Onufriev, A. V. (2016). Accuracy limit of rigid 3-point water
  models. *Journal of Chemical Physics*, 145(7), 074501.
  doi:10.1063/1.4960175. For `solvent.water.model: opc3` (`opc3.offxml`).
- Joung, I. S. and Cheatham, T. E. (2008). Determination of alkali and halide
  monovalent ion parameters for use in explicitly solvated biomolecular
  simulations. *Journal of Physical Chemistry B*, 112(30), 9020-9041.
  doi:10.1021/jp8001614. The Na+ and Cl- Lennard-Jones parameters of
  `openff-1.*` and `openff-2.*` are the TIP3P set of its Table 5. The
  openff-forcefields package keeps that table as `data/jc.csv`. The build
  cites it when it added Na+ or Cl- with one of these files.
- Sengupta, A., Li, Z., Song, L. F., Li, P. and Merz, K. M. (2021).
  Parameterization of monovalent ions for the OPC3, OPC, TIP3P-FB, and
  TIP4P-FB water models. *Journal of Chemical Information and Modeling*,
  61(2), 869-880. doi:10.1021/acs.jcim.0c01390. The Na+ and Cl- parameters
  of `opc3.offxml` (the Li/Merz 12-6 set for OPC3), which an OPC3 build uses
  instead of the solvent force field's ions.
- Jakalian, A., Bush, B. L., Jack, D. B. and Bayly, C. I. (2000). Fast,
  efficient generation of high-quality atomic charges. AM1-BCC model: I.
  Method. *Journal of Computational Chemistry*, 21(2), 132-146.
  doi:10.1002/(SICI)1096-987X(20000130)21:2<132::AID-JCC5>3.0.CO;2-P
- Jakalian, A., Jack, D. B. and Bayly, C. I. (2002). Fast, efficient
  generation of high-quality atomic charges. AM1-BCC model: II.
  Parameterization and validation. *Journal of Computational Chemistry*,
  23(16), 1623-1641. doi:10.1002/jcc.10128. Cite both AM1-BCC papers for
  `charges: {method: am1bcc}` and for a NAGL model of AM1-BCC charges.
- AshGC, the OpenFF NAGL model `openff-gnn-am1bcc-1.*.pt` (for example
  `openff-gnn-am1bcc-1.0.0.pt`). Cite the Sage 2.3.0 paper above, which
  describes it, and both AM1-BCC papers. "Methods to cite" names it
  `AshGC (OpenFF NAGL model <file>)`.
- Other OpenFF NAGL models,
  <https://github.com/openforcefield/openff-nagl-models>. Models before 1.0
  have no journal article; cite the model file the build names and its
  release.
- Wang, Y., Pulido, I., Takaba, K., Kaminow, B., Scheen, J., Wang, L. and
  Chodera, J. D. (2024). EspalomaCharge: machine learning-enabled ultrafast
  partial charge assignment. *Journal of Physical Chemistry A*, 128(20),
  4160-4167. doi:10.1021/acs.jpca.4c01287. For `charges: {method: espaloma}`.

## Libraries we build on

PolyzyMD does not reimplement trajectory input and output, neighbour searching,
secondary structure assignment or the standard statistical tests, so the
packages that provide them are cited as methods rather than as dependencies.

- Michaud-Agrawal, N., Denning, E. J., Woolf, T. B., and Beckstein, O. (2011).
  MDAnalysis: a toolkit for the analysis of molecular dynamics simulations.
  Journal of Computational Chemistry 32:2319-2327. doi:10.1002/jcc.21787
- Gowers, R. J., Linke, M., Barnoud, J., Reddy, T. J. E., Melo, M. N.,
  Seyler, S. L., Domanski, J., Dotson, D. L., Buchoux, S., Kenney, I. M., and
  Beckstein, O. (2016). MDAnalysis: a Python package for the rapid analysis of
  molecular dynamics simulations. Proceedings of the 15th Python in Science
  Conference, 98-105. doi:10.25080/Majora-629e541a-00e
- Smith, P., Ziolek, R. M., Gazzarrini, E., Owen, D. M., and Lorenz, C. D.
  (2019). On the interaction of hyaluronic acid with synovial fluid lipid
  membranes. Physical Chemistry Chemical Physics 21:9845-9857.
  doi:10.1039/C9CP01532A (we do not build on this paper. MDAnalysis asks that
  users of its HydrogenBondAnalysis cite it, and `functions.hydrogen_bonds` and
  `functions.hbond_count` use that class.)
- McGibbon, R. T., Beauchamp, K. A., Harrigan, M. P., Klein, C., Swails, J. M.,
  Hernandez, C. X., Schwantes, C. R., Wang, L.-P., Lane, T. J., and
  Pande, V. S. (2015). MDTraj: a modern open library for the analysis of
  molecular dynamics trajectories. Biophysical Journal 109:1528-1532.
  doi:10.1016/j.bpj.2015.08.015
- Virtanen, P., Gommers, R., Oliphant, T. E., and the SciPy 1.0 contributors
  (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python.
  Nature Methods 17:261-272. doi:10.1038/s41592-019-0686-2

## Structural algorithms

These define the geometry PolyzyMD measures, so the criteria and radii they set
decide what counts as a hydrogen bond, a helix or an exposed atom.

- Kabsch, W. (1976). A solution for the best rotation to relate two sets of
  vectors. Acta Crystallographica A32:922-923. doi:10.1107/S0567739476001873
- Theobald, D. L. (2005). Rapid calculation of RMSDs using a quaternion-based
  characteristic polynomial. Acta Crystallographica A61:478-480.
  doi:10.1107/S0108767305015266
- Liu, P., Agrafiotis, D. K., and Theobald, D. L. (2010). Fast determination of
  the optimal rotational matrix for macromolecular superpositions. Journal of
  Computational Chemistry 31:1561-1563. doi:10.1002/jcc.21439
- Kabsch, W. and Sander, C. (1983). Dictionary of protein secondary structure:
  pattern recognition of hydrogen-bonded and geometrical features. Biopolymers
  22:2577-2637. doi:10.1002/bip.360221211
- Shrake, A. and Rupley, J. A. (1973). Environment and exposure to solvent of
  protein atoms. Lysozyme and insulin. Journal of Molecular Biology 79:351-371.
  doi:10.1016/0022-2836(73)90011-9
- Bondi, A. (1964). van der Waals volumes and radii. Journal of Physical
  Chemistry 68:441-451. doi:10.1021/j100785a001
- Tien, M. Z., Meyer, A. G., Sydykova, D. K., Spielman, S. J., and Wilke, C. O.
  (2013). Maximum allowed solvent accessibilities of residues in proteins. PLoS
  ONE 8:e80635. doi:10.1371/journal.pone.0080635
- Arunan, E., Desiraju, G. R., Klein, R. A., Sadlej, J., Scheiner, S.,
  Alkorta, I., Clary, D. C., Crabtree, R. H., Dannenberg, J. J., Hobza, P.,
  Kjaergaard, H. G., Legon, A. C., Mennucci, B., and Nesbitt, D. J. (2011).
  Definition of the hydrogen bond (IUPAC Recommendations 2011). Pure and
  Applied Chemistry 83:1637-1641. doi:10.1351/PAC-REC-10-01-02
- Kuzmanic, A. and Zagrovic, B. (2010). Determination of ensemble-average
  pairwise root mean-square deviation from experimental B-factors. Biophysical
  Journal 98:861-871. doi:10.1016/j.bpj.2009.11.011 (not implemented)

## Sampling and uncertainty

Molecular dynamics frames are correlated, so the number of frames is not the
number of samples, and these sources say how to count the samples and how to
turn that count into an error bar.

- Grossfield, A., Patrone, P. N., Roe, D. R., Schultz, A. J., Siderius, D. W.,
  and Zuckerman, D. M. (2018). Best practices for quantifying sampling quality
  and uncertainty in molecular simulations. Living Journal of Computational
  Molecular Science 1:5067. doi:10.33011/livecoms.1.1.5067
- Hurlbert, S. H. (1984). Pseudoreplication and the design of ecological field
  experiments. Ecological Monographs 54:187-211. doi:10.2307/1942661 (the
  experimental unit: trends over a condition factor are fitted on condition
  means)
- Lazic, S. E. (2010). The problem of pseudoreplication in neuroscientific
  studies: is it affecting your analysis? BMC Neuroscience 11:5.
  doi:10.1186/1471-2202-11-5
- Grossfield, A. and Zuckerman, D. M. (2009). Quantifying uncertainty and
  sampling quality in biomolecular simulations. Annual Reports in Computational
  Chemistry 5:23-48. doi:10.1016/S1574-1400(09)00502-7
- Chodera, J. D., Swope, W. C., Pitera, J. W., Seok, C., and Dill, K. A.
  (2007). Use of the weighted histogram analysis method for the analysis of
  simulated and parallel tempering simulations. Journal of Chemical Theory and
  Computation 3:26-41. doi:10.1021/ct0502864
- Shirts, M. R. and Chodera, J. D. (2008). Statistically optimal analysis of
  samples from multiple equilibrium states. Journal of Chemical Physics
  129:124105. doi:10.1063/1.2978177
- Janke, W. (2002). Statistical analysis of simulations: data correlations and
  error estimation. In Quantum Simulations of Complex Many-Body Systems, NIC
  Series 10:423-445.
- Sokal, A. (1997). Monte Carlo methods in statistical mechanics: foundations
  and new algorithms. In Functional Integration, NATO ASI Series 361:131-192.
  doi:10.1007/978-1-4899-0319-8_6
- Lyman, E. and Zuckerman, D. M. (2007). On the structural convergence of
  biomolecular simulations by determination of the effective sample size.
  Journal of Physical Chemistry B 111:12876-12882. doi:10.1021/jp073061t
- Zhang, X., Bhatt, D., and Zuckerman, D. M. (2010). Automated sampling
  assessment for molecular simulations using the effective sample size. Journal
  of Chemical Theory and Computation 6:3048-3057. doi:10.1021/ct1002384
- Joint Committee for Guides in Metrology (2008). Evaluation of measurement
  data: guide to the expression of uncertainty in measurement. JCGM 100:2008.
- Flyvbjerg, H. and Petersen, H. G. (1989). Error estimates on averages of
  correlated data. Journal of Chemical Physics 91:461-466.
  doi:10.1063/1.457480 (block averaging, not implemented)
- Efron, B. (1979). Bootstrap methods: another look at the jackknife. Annals of
  Statistics 7:1-26. doi:10.1214/aos/1176344552 (not implemented)
- Efron, B. and Tibshirani, R. J. (1993). An introduction to the bootstrap.
  Chapman and Hall, New York. (not implemented)

## Inferential statistics

Conditions are compared across a handful of replicates with unequal variances,
which is what these tests and corrections are built for.

- Welch, B. L. (1947). The generalization of Student's problem when several
  different population variances are involved. Biometrika 34:28-35.
  doi:10.1093/biomet/34.1-2.28
- Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate:
  a practical and powerful approach to multiple testing. Journal of the Royal
  Statistical Society Series B 57:289-300.
  doi:10.1111/j.2517-6161.1995.tb02031.x
- Cohen, J. (1988). Statistical power analysis for the behavioral sciences,
  2nd edition. Lawrence Erlbaum Associates, Hillsdale.
- Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect
  size and related estimators. Journal of Educational Statistics 6:107-128.
  doi:10.3102/10769986006002107 (the small-sample correction behind `hedges_g`)
- Scott, D. W. (1992). Multivariate density estimation: theory, practice, and
  visualization. Wiley, New York. doi:10.1002/9780470316849

## Convergence

A trajectory that has not relaxed gives a confident wrong answer, so these
sources define where production starts and whether the runs agree.

- Chodera, J. D. (2016). A simple method for automated equilibration detection
  in molecular simulations. Journal of Chemical Theory and Computation
  12:1799-1805. doi:10.1021/acs.jctc.5b00784 (not implemented)
- Yang, W., Bitetti-Putzer, R., and Karplus, M. (2004). Free energy
  simulations: use of reverse cumulative averaging to determine the
  equilibrated region and the time required for convergence. Journal of
  Chemical Physics 120:2618-2628. doi:10.1063/1.1638996 (not implemented)
- Klimovich, P. V., Shirts, M. R., and Mobley, D. L. (2015). Guidelines for the
  analysis of free energy calculations. Journal of Computer-Aided Molecular
  Design 29:397-411. doi:10.1007/s10822-015-9840-9
- Hess, B. (2002). Convergence of sampling in protein simulations. Physical
  Review E 65:031910. doi:10.1103/PhysRevE.65.031910
- Knapp, B., Ospina, L., and Deane, C. M. (2018). Avoiding false positive
  conclusions in molecular simulation: the importance of replicas. Journal of
  Chemical Theory and Computation 14:6127-6138. doi:10.1021/acs.jctc.8b00391
