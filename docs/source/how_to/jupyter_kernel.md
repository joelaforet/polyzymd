# Use a PolyzyMD environment in Jupyter or VS Code

The notebooks that `polyzymd project init` puts in `notebooks/` need a
PolyzyMD environment to run. This page shows three ways to give them one:
start JupyterLab from the environment, register the environment as a Jupyter
kernel, or pick its Python in VS Code.

The commands use the `build` environment, which has JupyterLab and the
analysis tools. The `analysis` environment works the same way. Change
`~/polyzymd` to the folder of your PolyzyMD clone.

## Start JupyterLab from the environment

The simplest way needs no setup. Start JupyterLab with `pixi run`, from the
folder of the notebook:

```bash
cd ~/pz_quickstart/notebooks
pixi run --manifest-path ~/polyzymd/pixi.toml -e build jupyter lab --no-browser
```

Open the link with a token that JupyterLab prints. Every notebook it opens runs
in the `build` environment.

## Register the environment as a kernel

Register it once if you start Jupyter in another way, such as a JupyterLab
that you installed yourself or a JupyterHub on a cluster:

```bash
pixi run --manifest-path ~/polyzymd/pixi.toml -e build \
  python -m ipykernel install --user --name polyzymd-build --display-name "PolyzyMD (build)"
```

The output names the folder of the new kernel:

```text
Installed kernelspec polyzymd-build in /home/me/.local/share/jupyter/kernels/polyzymd-build
```

Then choose **PolyzyMD (build)** in the kernel menu of a notebook. Check the
list of kernels with `jupyter kernelspec list`.

The kernel starts the environment's Python directly, without activating the
environment, so the programs of the environment are not on its `PATH`.
Analyses and the notebooks of `notebooks/` do not need them. To build a system
in a notebook, which runs PACKMOL, start JupyterLab from the environment as
above.

The kernel points at the environment's folder in the clone. After you move or
delete the clone, register the kernel again. To remove it:

```bash
jupyter kernelspec remove polyzymd-build
```

## Use the environment in VS Code

VS Code runs a notebook with any Python interpreter. Open the notebook, click
**Select Kernel** in the top right, choose **Python Environments**, and pick
the environment's Python:

```text
~/polyzymd/.pixi/envs/build/bin/python
```

If it is not listed, choose **Select Another Kernel**, then
**Python Environments**, then **Enter interpreter path**, and type the path.
On Windows with WSL, open the folder with the WSL extension of VS Code first,
so that VS Code runs on the Linux side, where the environment is.

## Related

- {doc}`own_analysis`: what `load_replicates.ipynb` does.
- {doc}`../tutorials/star_sited_monomers`: the notebook on monomers with `*` sites.
