# Set Up a Contributor Environment

Use this page when you want to modify PolyzyMD, run tests, build the docs, or
work on new analysis plugins.

## What This Environment Is For

The contributor setup gives you:

- an editable PolyzyMD install through `pixi.toml`
- the test and docs toolchain
- the same dependency story used by the project itself

<!-- IMAGE OPPORTUNITY: Add a short terminal transcript showing `pytest`,
`ruff check`, and `make -C docs clean html` succeeding inside the build
environment. -->

## 1. Clone the Repository

```bash
git clone https://github.com/joelaforet/polyzymd.git
cd polyzymd
```

## 2. Install the Build Environment

```bash
pixi install -e build
```

PolyzyMD is installed in editable mode by default, so source changes are picked
up immediately. If you want an interactive shell inside the environment, run
`pixi shell -e build`.

## 3. Run the Core Contributor Checks

```bash
pixi run -e build pytest tests/ -x -q --tb=short --ignore=tests/utils/test_packmol.py
pixi run -e build ruff check src/
pixi run -e build black --check src/
pixi run -e build make -C docs clean html
```

## 4. Optional Tooling Notes

### OpenEye Toolkit

If you have an OpenEye license and want the commercial charging backend:

```bash
pixi run -e build python -m pip install openeye-toolkits
```

If you prefer an interactive shell, enter it first with `pixi shell -e build`
and then run `python -m pip install openeye-toolkits` inside that shell. Avoid
installing OpenEye into a system Python that PolyzyMD will not use.

### CUDA Environments for Cluster Work

If you are debugging or developing on GPU clusters, install the CUDA-specific
pixi environment that matches the target cluster driver version:

```bash
pixi install -e cuda-12-6
pixi shell -e cuda-12-6
```

## 5. Keep Your Environment Updated

```bash
git pull
pixi install -e build
```

## Related Pages

- {doc}`index`
- {doc}`packaging`
- {doc}`contributing`
- {doc}`extending_analyses`
