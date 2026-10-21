"""Building blocks the analyses are made of: loading, windows, statistics and plotting.

Sub-modules
-----------
loader
    Trajectory loading, time parsing, frame conversion.
window
    The production window of a replicate after the equilibration time.
selections
    Extended selection syntax (midpoint, COM), position retrieval.
diagnostics
    Selection diagnostics, equilibration validation.
centroid
    The frame closest to the iterative average structure, for centroid references.
topology
    Checks that a topology carries the bonds an analysis needs.
statistics
    Mean, standard error and Student t interval of replicate values.
inferential_statistics
    t tests, effect sizes and the Benjamini-Hochberg correction.
autocorrelation
    Statistical inefficiency, effective sample size and detected equilibration.
aa_classification
    Maximum accessible surface area of each amino acid.
groupings
    Physicochemical classes of amino acids.
plotting
    Axis styling, condition colours, grouped bars, uncertainty bands, figure saving.
"""
