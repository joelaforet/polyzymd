# Analysis APIs

Use these pages for the analysis plugin framework, its base classes, and the
shared analysis utilities. No shipped analysis is a plugin any more: every
analysis that `polyzymd analyze` offers is a function in
`polyzymd.analyses.functions` run through the study API, described in
{doc}`../explanation/analysis_api` and {doc}`../reference/analysis_functions`.
The plugin framework is being removed.

```{toctree}
:maxdepth: 1

Analysis Plugin Framework API <analyses>
Analysis Base Classes <analyses_base>
Analysis Shared Utilities <analyses_shared>
```
