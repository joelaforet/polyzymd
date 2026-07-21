# Development Workflow

## Project Priorities

Optimize PolyzyMD work in this order:

1. Scientific correctness and reproducibility.
2. Low-friction user experience.
3. OpenMM and GROMACS parity for engine-facing features.
4. Small, reviewable, maintainable changes.
5. Documentation that matches the code.

Joe retains authority over scientific assumptions and release readiness. Pause
and ask before choosing chemistry, topology, parameterization, force-field, or
scientific-interpretation assumptions. Do not silently substitute agent
judgment in these areas.

## Collaboration and Progress

- Use one writer per worktree. Give review, architecture, documentation, and
  validation agents read-only roles unless they receive an explicit writing
  handoff.
- Use separate worktrees for parallel writing tasks.
- Report milestone updates rather than every command. Useful milestones include
  scope confirmation, implementation completion, focused tests, engine-specific
  artifacts, parity checks, short dynamics, documentation, and independent
  review.
- Keep work steerable. When evidence changes the likely implementation path,
  report it before expanding the scope.

## Atomic Scope

Start from one end-user goal and define explicit non-goals. A pull request must
be a small, self-contained, independent change that addresses one well-defined
goal. Multiple atomic PRs may compose a larger feature.

- Below roughly 300 changed lines: normally proceed.
- Above 300 changed lines: inspect for scope creep, generated content, and
  excessive tests; split when independent behavior exists.
- Above 500 changed lines: pause and explain why the change cannot be divided.

Test code counts toward the review burden. Prefer the smallest test set that
proves the behavior and protects the scientific contract. Do not add layers of
repetitive fixtures, mocks, or abstractions without a concrete failure mode.

Treat `and` in a commit or PR subject as a signal to check whether two units of
work should be separated. It is a scope heuristic, not a grammatical ban.

## Validation and Engine Parity

Use a proportional validation ladder: static checks, focused tests, subsystem
tests, representative integration, scientific acceptance, broad regression,
and documentation. Expensive validation must follow cheaper gates.

Engine-facing features require OpenMM and GROMACS coverage. Match molecular
identity and intended mechanics: particles, masses, bonds, constraints,
charges, parameters, exclusions/pairs, coordinates, and box. Allow justified
numerical differences caused by engine algorithms, including small PME energy
differences; document the tolerance and evidence.

On `jlaforetlaptop`, query OpenMM's API for available platforms before choosing
CUDA. Do not assume that a platform is available from the machine name alone.

When visual molecular inspection remains necessary, hand Joe absolute paths to
the topology and, when available, trajectory. State what automation checked and
what still requires his PyMOL judgment.

## Commits and Publication

After an atomic change passes its required gates:

1. Review the complete owned diff.
2. Verify `git config user.name` and `git config user.email` identify Joe.
3. Create an atomic conventional commit using the configured identity.
4. Explain non-obvious decisions and scientific intent in the commit body.
5. Do not add Codex branding, generated-by text, or co-author trailers.
6. Push validated commits to the working branch automatically for authorized
   development tasks.

Do not commit or push incomplete, failing, or scientifically unresolved work.
Read-only inspection does not authorize commits or pushes. Never force-push a
published/shared branch without Joe's explicit approval; use
`--force-with-lease` only after that approval.

Do not open a pull request until Joe confirms development is complete. After
confirmation, Codex opens the PR and reports the user goal, rationale,
scientific assumptions, validation, documentation, artifact paths, limitations,
linked issues, and dependencies. Codex never merges PRs; Joe reads and merges
every PR manually.

The GitHub connector lacked PR-write permission, so I used the authenticated gh fallback. I will not merge it.

## Maintain Living Guidance

At the end of each change, compare every loaded repository instruction, Skill,
and reference against the behavior and tool access observed during the task.
When code, scientific contracts, commands, branch names, release policy,
documentation structure, tools, permissions, or recurring workflows change:

1. Update the affected repository guidance and personal PolyzyMD Skill in the
   same atomic task when access and scope permit.
2. Remove or revise obsolete branch aliases, commands, limitations, and
   workarounds instead of accumulating contradictory history.
3. Validate edited Skills with the Skill validator and run the relevant
   repository documentation or static gates.
4. If the guidance cannot be updated, report the exact stale file and proposed
   replacement text to Joe before closing the task.

Only persist evidence-backed, reusable knowledge. Do not turn a one-off failure
or system-specific tuning choice into a universal rule.

## Release Flow

- `main` contains released code.
- `release/1.3` stabilizes the v1.3 analysis release. Until migration, its
  current branch name is `feature/v1.3.0-rc5`.
- `release/1.4` integrates conjugation work. Until migration, its current branch
  name is `conjugation-engine-refactor`.
- Branch bug fixes from the current release line and merge them back through
  atomic PRs.
- Tag immutable beta snapshots as `vX.Y.Z-rc.N`.
- Forward-integrate release-line fixes into the next release line regularly.
- Run lightweight lint, import, unit, and documentation checks on ordinary PRs.
  Reserve GPU simulations and resource-intensive scientific regression for
  feature acceptance and major release gates.
