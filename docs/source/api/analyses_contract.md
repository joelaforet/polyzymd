# The observable contract

API reference for `polyzymd.analyses.contract`, the module a plugin author
imports. It defines what a plugin provides, the four observable kinds, how a
replicate reduces, how conditions are compared, and the factory that turns a
plugin into a runnable analysis.

`Observable` is one measured quantity from one replicate. `iter_frames` walks
the production window. `contract_analysis` builds the
`polyzymd.analyses.base.Analysis` subclass that discovery finds, and is the last
line of every plugin module.

`ObservableEstimate`, `ObservableAggregate` and `ObservableComparison` are what
the framework writes: one estimate per replicate, one aggregate per condition,
one comparison per tested observable per condition pair.

For the contributor walkthrough see
{doc}`../contributor_guide/analysis_plugins/index`.

```{eval-rst}
.. automodule:: polyzymd.analyses.contract
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```
