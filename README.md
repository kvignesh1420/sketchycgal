# sketchycgal
Implementation of the SketchyCGAL algorithm from https://arxiv.org/pdf/1912.02949.pdf

### Instructions

```bash
# download the datasets
$ bash prepare_data.sh

# build the C++ library and bind to python as well
$ bash build.sh

# C++ runs
$ bazel-bin/solver/solver <path to graph .mtx file>

# clean bazel cache
$ bazel clean --expunge
```

For python based experiments, please create a virtual environment

```bash
$ python3.9 -m virtualenv .venv
$ source .venv/bin/activate
$ pip install -r requirements.txt
# python runs (by default run on G67 graphs.)
# Change the path in code for other larger graphs.
(.venv) $ python experiments.py
```

### Profiling

Ensure that you have `perf` on your system

```bash

# an example for delauney20 graph
# record the data
$ OMP_NUM_THREADS=8 perf record bazel-bin/solver/solver data/DIMACS10/delaunay_n20/delaunay_n20/delaunay_n20.mtx

# report the data
$ perf report -i perf.data

```

### Singularity based execution on Greene

Please use the `cpu-hpc.sbatch` file to launch a jupyter slurm job. This will provide you with a link to jupyter environment, through which you can interactively build the library and explore the functionality. Please feel free to modify the resource allocations as per experimental requirements.

```bash
$ sbatch  cpu-hpc.sbatch
```