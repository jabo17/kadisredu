# KaDisRedu: Distributed Reductions for the Maximum Weight Independent Set Problem

This project is joint work between Jannick Borowitz[^1], Ernestine Großmann[^2], and Matthias Schimek[^3].

[![DOI](https://zenodo.org/badge/807071373.svg)](https://doi.org/10.5281/zenodo.17174407)

---

In the following, we describe the necessary steps to build and use our algorithms.
If you use our algorithms in an academic work, please acknowledge our work by citing the corresponding paper[^4].
If you are particularly interested in reproducing artifacts in the corresponding paper [^4] of this repository, we refer you our [reproducibility repository](https://github.com/jabo17/kadisredu-reproducibility).

## Dependencies
We tested the following software stack to compile this code:
- CMake version >= 3.25
- GCC 13
- IntelMPI 2021.11
- Intel TBB 2021.4 (used libraries also require Intel TBB)
- Google Sparse Hash (libsparsehash) (required by external libraries)

## Build our algorithms
Note that you require an internet connection because some external libraries are downloaded from the internet via `FetchContent` (CMake).

```bash
cd distributed-reductions
cmake --preset=Release -S . -B build
cmake --build build --target kadisredu_app --parallel 
```

## Run our algorithms

```bash
# your settings
graph="snap_soc-pokec-relationships-uniform" # from `examples/instances`
json_output="output.json"
algo="RGA.toml"
time_limit=7200 # seconds
cores=1

cd distributed-reductions
I_MPI_PIN=1 I_MPI_PIN_DOMAIN=core I_MPI_PIN_ORDER=compact I_MPI_JOB_TIMEOUT=$timelimit mpiexec -n $cores ./build/apps/kadisredu_app --time_limit ${time_limit} --seed 0 --warmup_mpi --json_output_path "${json_output}" --kagen_option_string "file;filename=../examples/instances/${graph}.parhip;distribution=balance-edges" --configs "config/${algo}"
```

### Output
We output the measured running times beside other metrics to a json file `$json_output`.
The following fields might be of interest:
- overall running time: `.timer.KaDisRedu.statistics.max[0]`
- overall solution weight: `.solver.solution_weight`
- reduce time: `.timer.KaDisRedu.reduce.statistics.max[0]`

### Algorithms
We provide for each solver of the paper a toml configuration in `distributed-reductions/configs`.
If you are only interested in our reduce algorithms `DisReduA` or `DisReduS`, you can use `RGA` or `RGS`, respecitvely.
To partition the graph first, before reducing it, use the configurations `PRGA` or `PRGS`.
In the output you can find the metrics of the reduce phase (`.solver.R0`).
You can print the reduced graph in a file in the METIS format by adding `--kernel <file>` to the command.
Similarly, you can add `--independent_set <file>` to print the maximal inpendent set.
The i-th line indicates whether the i-th vertex is in the solution (1) or not (0).

## References
[^1]: Jannick Borowitz: [jannick.borowitz@kit.edu](mailto:jannick.borowitz@kit.edu) (Karlsruhe Institute of Technology)
[^2]: Ernestine Großmann: [e.grossmann@informatik.uni-heidelberg.de](mailto:e.grossmann@kit.uni-heidelberg.de) (Heidelberg University)
[^3]: Matthias Schimek: [schimek@kit.edu](mailto:schimek@kit.edu) (Karlsruhe Institute of Technology)
[^4]: TODO Link to paper
