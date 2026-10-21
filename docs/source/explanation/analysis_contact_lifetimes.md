# How long polymer contacts last

This page explains what `polyzymd analyze contacts --run mean_lifetime`
measures, why it uses the Kaplan-Meier estimator rather than a plain mean, and
which choices change the number.

## Contact events

A residue is in contact with the polymer on a frame by one of the two
definitions of {doc}`../how_to/analysis_contacts_quickstart`: the polymer
buries it (`method=occlusion`) or comes within a cutoff of it
(`method=distance`). An **event** is a run of consecutive production frames in
which one residue is in contact. A run of `k` frames lasts `k` times the
spacing of the frames. The events of every measured residue are pooled, so
the result answers: once the polymer touches a residue, how long does the
contact typically last? For a monomer type, `<type>_mean_lifetime`, a contact
is one with that type's atoms alone.

Events are found per residue, whichever polymer molecule or monomer makes the
contact, so a residue passed from one chain to another without a break stays
one event.

## Censoring and the Kaplan-Meier estimator

An event still under way at the last production frame, or already under way
at the first, has an unknown length: it lasted at least as long as observed.
Such an event is **censored**. Long events are the most likely to be cut off,
so either counting censored events as finished or leaving them out makes a
plain mean too short.

The Kaplan-Meier estimator (Kaplan and Meier 1958) is the standard way to
estimate a distribution of durations when some durations are only lower
bounds, as for patients still alive when a clinical study ends. It estimates
the survival function S(t), the probability that an event lasts longer than
`t`, using every event, finished or censored. PolyzyMD computes it with
`scipy.stats.ecdf` on `scipy.stats.CensoredData`.

The reported `mean_lifetime` is the **restricted mean** (Royston and Parmar
2013): the area under S(t) from 0 to the time the production frames span. It
is the mean duration with every duration cut at that span, so it is defined
even when most events outlast the trajectory, and it is comparable between
replicates of the same length. When the longest event is censored, S(t) stays
at its last value up to the span, as the Kaplan-Meier estimate does beyond the
last observation. `censored_fraction` reports how many events were censored;
when it is large, the restricted mean mostly reflects the trajectory length.
`lifetime_events` reports how many events there were.

Each replicate gives one restricted mean, and conditions are compared with the
replicate as the sampling unit, as for every other result.

## Choices that change the number

**Frame spacing.** A contact shorter than the spacing of the frames is
invisible, and a break shorter than the spacing is not seen. Measuring every
tenth frame therefore gives fewer, longer events. On a 363 K lipase A
replicate with a 50:50 SBMA-EGMA polymer and distance contacts, frames 40 ps
apart gave 6,961 events and a restricted mean of 0.36 ns. Compare conditions
at the same `--stride` and frame spacing only.

**Tolerance for brief breaks.** With `--set tolerance_ps=...`, absences of at
most that many picoseconds between two contacts are filled, with MDAnalysis
`lib.correlations.correct_intermittency` (Gowers and Carbone 2015), so they do
not end the event. The tolerance is converted to whole frames of each
replicate's spacing. On the same replicate, a tolerance of 40 ps (one frame)
raised the restricted mean from 0.36 to 0.86 ns, and 120 ps to 1.87 ns. Laage
and Hynes (2008) showed that residence times depend strongly on such a
tolerance, so the default is 0, and a tolerance, when used, should be reported
with the result.

**Contact definition.** Occlusion and distance contacts, and their thresholds
and cutoffs, count different frames as contacts and so give different events.

## References

**Kaplan EL, Meier P.** (1958) "Nonparametric estimation from incomplete
observations." *J Am Stat Assoc* 53:457-481.
https://doi.org/10.1080/01621459.1958.10501452

**Royston P, Parmar MKB.** (2013) "Restricted mean survival time: an
alternative to the hazard ratio for the design and analysis of randomized
trials with a time-to-event outcome." *BMC Med Res Methodol* 13:152.
https://doi.org/10.1186/1471-2288-13-152

**Laage D, Hynes JT.** (2008) "On the residence time for water in a solute
hydration shell: application to aqueous halide solutions." *J Phys Chem B*
112:7697-7701. https://doi.org/10.1021/jp802033r

**Gowers RJ, Carbone P.** (2015) "A multiscale approach to model hydrogen
bonding: the case of polyamide." *J Chem Phys* 142:224907.
https://doi.org/10.1063/1.4922445

## See Also

- [Contacts Quick Start](../how_to/analysis_contacts_quickstart.md) — commands and settings
- [Shipped analysis functions](../reference/analysis_functions.md) — what each function measures
