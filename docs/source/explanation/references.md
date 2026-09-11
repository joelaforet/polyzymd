# Methods and references

Every method PolyzyMD's analysis module implements or adapts comes from a
published source, and this page lists them. If you add a method, add its
reference here and in a NumPy-style `References` section in the module
docstring. Entries marked as not implemented describe planned work and are
listed so nobody cites them from code that does not yet run them.

To cite PolyzyMD itself, see `CITATION.cff` in the repository root.

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
  users of its HydrogenBondAnalysis cite it, and the hydrogen bonds plugin uses
  that class.)
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
- Tukey, J. W. (1949). Comparing individual means in the analysis of variance.
  Biometrics 5:99-114. doi:10.2307/3001913
- Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate:
  a practical and powerful approach to multiple testing. Journal of the Royal
  Statistical Society Series B 57:289-300.
  doi:10.1111/j.2517-6161.1995.tb02031.x
- Cohen, J. (1988). Statistical power analysis for the behavioral sciences,
  2nd edition. Lawrence Erlbaum Associates, Hillsdale.
- Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect
  size and related estimators. Journal of Educational Statistics 6:107-128.
  doi:10.3102/10769986006002107 (small-sample correction, not implemented)
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
