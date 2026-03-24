# Set Up a Contributor Environment

Use this page when you want to modify PolyzyMD, run tests, build the docs, or
work on new analysis and comparison features.

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

## 2. Install and Activate the Build Environment

```bash
pixi install -e build
pixi shell -e build
```

PolyzyMD is installed in editable mode by default, so source changes are picked
up immediately.

## 3. Run the Core Contributor Checks

```bash
pytest tests/ -x -q --tb=short --ignore=tests/test_packmol.py
ruff check src/
black --check src/
make -C docs clean html
```

## 4. Optional Tooling Notes

### OpenEye Toolkit

If you have an OpenEye license and want the commercial charging backend:

```bash
pixi shell -e build
pip install openeye-toolkits
```

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
- {doc}`../tutorials/packaging`
- {doc}`../tutorials/contributing`
- {doc}`../tutorials/extending_comparators`
- {doc}`../tutorials/extending_plotters`
