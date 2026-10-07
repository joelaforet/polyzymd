# How the study API analyzes a set of simulations

The study API runs your analysis code on every replicate of every condition
of a study. It stores each result with a record of its inputs, and compares
the conditions with the replicate as the sampling unit. This page explains
what the API does and why. To use it, see {doc}`../how_to/study_api`. For
every signature, see {doc}`../reference/study_api`.

The shipped analyses, such as the radius of gyration or the hydrogen-bond
occupancy, are ordinary functions that use the same interface. So your own
code can do anything that they do.

## What you write and what PolyzyMD does

| You write | PolyzyMD does |
|---|---|
| A function that measures one frame, or one replicate | Loads each replicate as a `Universe`, with its production segments in order and its equilibration window removed |
| The atom selections of the function, as selection strings | Builds the `AtomGroup` of each selection in each replicate, and runs the function on the production frames |
| Optionally, how the frames of one replicate become one value | Stores each per-frame and per-replicate result with the code, arguments and input files that produced it, and reuses it only when all of them match |
| Nothing more for the statistics | Computes the mean, the standard error and the 95 % Student t interval across replicates, tests each condition against the control, and corrects for multiple tests |

The per-frame interface follows MDAnalysis
[`AnalysisFromFunction`](https://docs.mdanalysis.org/stable/documentation_pages/analysis/base.html).
Your function receives `AtomGroup`s at the current frame and returns a
number. An `AtomGroup` belongs to one `Universe`. So you give a selection
string, `pz.select("protein")`, and PolyzyMD builds the `AtomGroup` again for
each replicate. It records the string.

## Why replicates line up by label

A function can return one value per label, such as one value per residue.
PolyzyMD aggregates the values by label, never by position. Residue 57 of
one replicate is always compared with residue 57 of another, even when a
replicate has a different number of entries. A label that one replicate
lacks is an error, unless you say which value it takes. This prevents a
silent shift of a whole profile by one residue.

## Why references and files are recorded by content

PolyzyMD builds a reference structure once per replicate, in a separate
universe. So the trajectory never changes.

The record of a result names each file that an argument points to by its
name and its SHA-256, not by its location. So an edit to a reference file
measures the replicates again, and a move of the study folder does not.

## The statistics

The statistics follow Grossfield et al. (2018), the LiveCoMS best practices
for the uncertainty of molecular simulations.

**The replicate is the sampling unit.** The frames of one replicate are
correlated, so their number never makes an interval narrower. PolyzyMD
reports the {term}`statistical inefficiency` g and the {term}`n_eff` of each
replicate's time series, computed with `pymbar.timeseries`, as diagnostics.

**One equilibration window for all replicates.** The window that you set
applies to every replicate of every condition. For each time series, PolyzyMD
also reports where `pymbar.timeseries.detect_equilibration` (Chodera 2016)
finds the start of the equilibrated region. If that point is after your
window, a warning says that the window can be too short for that replicate.
PolyzyMD never changes the window. See {doc}`convergence_detection`.

**The interval.** The 95 % confidence interval is the mean plus or minus the
Student t factor for n replicates times the standard error.

**The test.** PolyzyMD uses Welch's t test by default, or Student's t test,
with the {term}`Benjamini-Hochberg` correction.

**The correction family is one outcome.** For one number per replicate, the
family is the conditions compared with the control for one quantity, such as
the mean distance of one pair. A family is the set of tests behind one
conclusion (Bender and Lange 2001; Rubin 2021). A false discovery rate that is
controlled in separate families stays controlled overall (Benjamini and
Yekutieli 2001). A pool of unrelated outcomes can hide real effects or
inflate weak ones (Efron 2008).

For a labelled result, such as a per-residue profile, the family is every
label of every compared condition. A scan of a profile for the labels that
changed is a search over that whole set. A false discovery rate holds for the
discoveries of a search only when the family is the set searched (Benjamini
2010). The text report gives, for each condition, the number of labels that
are significantly lower and higher, and lists them. The JSON report keeps
every row. Each comparison line gives the raw p value, the adjusted p value
and the family size, `family <m>`.

If you make one claim from several outcomes together, such as "any of these
pairs changed", put those outcomes in one family. Pass their p values to
`polyzymd.analyses.shared.inferential_statistics.benjamini_hochberg`, and
report that family.

**Untestable rows.** A row is not testable when a condition has fewer than
two replicates, or when both conditions have the same value in every
replicate, such as a residue that is never in contact. Such a row takes no
part in the correction. You can also leave out a label whose value the other
labels fix, such as the last of a set of fractions that sum to 1.

**Bounded quantities.** A quantity with a strict lower or upper limit is not
Gaussian, and a t interval can extend past the limit. Grossfield et al.
recommend the bootstrap for such quantities, but also say that bootstrap
intervals are unreliable for small samples. Three to five replicates is the
usual case. So PolyzyMD keeps the t interval, and warns when it extends past
the `bounds` of the result. When every replicate has the same value,
PolyzyMD reports the interval as not estimable, not as an interval of zero
width.

## The figures

Grossfield et al. recommend that a figure shows every point when there are
fewer than 10 independent measurements. So every PolyzyMD figure shows each
replicate value, beside the mean and the interval.

Each figure with error bars or a band has a footnote. The footnote says what
the bars are, of which quantity, and across which replicates. An example is
"Error bars: 95% Student t confidence interval of the condition mean across
n = 5 replicates; production window t >= 10ns." When the conditions have
different numbers of replicates, the footnote says that n is given per
condition.

A distribution figure draws a Gaussian kernel density estimate of the frame
values. Near a finite bound of the series, PolyzyMD corrects the estimate by
reflection (Schuster 1985; Silverman 1986, section 2.10). So the curve stays
within the physical range and still integrates to 1.

Figures read the stored results, so they read no trajectory.

## Why stored results are reused

A per-replicate result is stored with a record of the function, the
arguments, the config hash, the input files, the window, the frames and the
software versions. PolyzyMD reuses the result only when every field except
the versions and the bounds matches the new call. A change to the function,
an argument or an input file, or a longer trajectory, measures the affected
replicates again.

A function defined in a notebook cell, or a lambda, has no stable source
file. PolyzyMD then hashes its bytecode, and warns that the record
identifies the function less reliably.

## What is out of scope

- An analysis that needs all replicates at once inside the measurement, such
  as a principal component basis, clustering or a Markov state model. Load
  the universes with `replicate.universe()` and fit such models yourself.
- A quantity defined between two conditions, not for one replicate, such as
  a divergence between two distributions. Compute it in your own script from
  the stored values (`replicate_table`) or from the universes.
- Intervals that are exact for a non-Gaussian quantity from three to five
  replicates. No method gives these. So PolyzyMD warns when an interval
  extends past a limit.

## Shipped analyses

The shipped analyses are functions in `polyzymd.analyses.functions` that use
this interface. `polyzymd analyze NAME` runs them from the command line. Each
one calls MDAnalysis or MDTraj for the per-frame measurement. For the list,
with the measurement and the reduction of each one, see
{doc}`../reference/analysis_functions`.

## References

Bender, R.; Lange, S. Adjusting for Multiple Testing: When and How? *J. Clin.
Epidemiol.* **2001**, 54 (4), 343-349. doi:10.1016/S0895-4356(00)00314-0

Benjamini, Y. Discovering the False Discovery Rate. *J. R. Stat. Soc. B*
**2010**, 72 (4), 405-416. doi:10.1111/j.1467-9868.2010.00746.x

Benjamini, Y.; Yekutieli, D. The Control of the False Discovery Rate in
Multiple Testing under Dependency. *Ann. Stat.* **2001**, 29 (4), 1165-1188.
doi:10.1214/aos/1013699998

Chodera, J. D. A Simple Method for Automated Equilibration Detection in
Molecular Simulations. *J. Chem. Theory Comput.* **2016**, 12 (4), 1799-1805.
doi:10.1021/acs.jctc.5b00784

Efron, B. Simultaneous Inference: When Should Hypothesis Testing Problems Be
Combined? *Ann. Appl. Stat.* **2008**, 2 (1), 197-223. doi:10.1214/07-AOAS141

Grossfield, A.; Patrone, P. N.; Roe, D. R.; Schultz, A. J.; Siderius, D. W.;
Zuckerman, D. M. Best Practices for Quantifying Uncertainty and Sampling
Quality in Molecular Simulations. *Living J. Comput. Mol. Sci.* **2018**, 1 (1),
5067. doi:10.33011/livecoms.1.1.5067

Rubin, M. When to Adjust Alpha during Multiple Testing: A Consideration of
Disjunction, Conjunction, and Individual Testing. *Synthese* **2021**, 199,
10969-11000. doi:10.1007/s11229-021-03276-4

Schuster, E. F. Incorporating Support Constraints into Nonparametric Estimators
of Densities. *Commun. Stat. Theory Methods* **1985**, 14 (5), 1123-1136.
doi:10.1080/03610928508828965

Silverman, B. W. *Density Estimation for Statistics and Data Analysis*; Chapman
and Hall: London, 1986; Section 2.10.
