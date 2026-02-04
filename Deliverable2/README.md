# PARCO-Computing-2026-242734 
# Deliverable 2

## Author
[Sandrini Leonardo] (https://github.com/Suge42)

Distributed sparse matrix-vector multiplication (SpMV) benchmark using MPI parallelization.

This project implements a distributed-memory parallel SpMV, with **Row-wise matrix distribution** across MPI ranks using **MPI Inter-node Communication**.

## Project Structure

```
Deliverable2/
├── src/                            # C source code
│   ├── execute_mpi_generating.c    # Main MPI SpMV program for weak scaling testing
│   ├── execute_mpi_reading.c       # Main MPI SpMV program for strong scaling testing
│   ├── libraries/                      # Additional C code used
│       ├── SpMV.c                      # Function to perform and check SpMV
│       ├── data_management.c           # General function for data collection
│       ├── bubblesort.c                # Bubblesort function for COO to CSR convertion
│       ├── generator.c/h               # Generator functions for weak scaling
│       ├── matrix_reading.c/h          # Matrix reading and conversion functions for strong scaling
│       ├── mmio.c                      # Library for matrix market reading
│       └── *.h                         # Header files for previous .c
│
├── matrices/                       # Sparse matrices (.mtx) - download required
│   └── *.mtx                       # Various matrices used in testing
│
├── scripts/                        # matrix download, PBSs and plotting scripts
│   ├── matrix_download.sh          # Download desired matrix from Suitsparse
│   ├── del2_weak_scaling.pbs       # PBS script for strong scaling testing
│   ├── load_balance_sweep.sh       # PBS script for manual testing of "execute_mpi_generating.c"
│   ├── del2_strong_scaling.pbs     # PBS script for strong scaling testing
│   ├── del2_reading.pbs            # PBS script for manual testing of "execute_mpi_reading.c"
│   ├── plotter_ws.py               # Python script for weak scaling data plotting
│   └── plotter_ss.py               # Python script for strong scaling data plotting
│
├── results/                        # Benchmark output data
│   ├── *.err                       # Error result files
│   ├── *.out                       # Main output files, with workflow and data
│   └── to_plot/                        # Results to be plotted
│       └── *.txt                       # Result files showing only the data, without the workflow
│
├── plots/                          # Generated figures (PNG)
│   ├── del2_ws_speedup_vs_processes.png
│   ├── del2_ws_speedup_vs_matrix_size.png
│   └── del2_ss_speedup_vs_processes.png
│
└── README.md                       # This file (project documentation)
```

## Quick Start

### Prerequisites
- **Compiler**: GCC ≥ 9.1
- **System**: Linux cluster with multi-node capability
- **Memory**: ≥ 32GB RAM per node recommended
- **MPI**: MPICH 3.2.1 or compatible MPI implementation
- **Python** (optional): For analysis scripts (matplotlib, pandas, numpy)

### Manual Compilation

Important to enclude each useful *.c file in compiling, since header files are used in code

```bash
# Compile matrix generating executable
mpicc -O2 -g -Wall -Wextra \
  ./src/execute_mpi_generating.c \
  ./src/libraries/SpMV.c \
  ./src/libraries/data_management.c \
  ./src/libraries/bubblesort.c \
  ./src/libraries/generator.c \
  -o del2_g

# Compile matrix reading executable
mpicc -O2 -g -Wall -Wextra \
  ./src/execute_mpi_reading.c \
  ./src/libraries/SpMV.c \
  ./src/libraries/data_management.c \
  ./src/libraries/bubblesort.c \
  ./src/libraries/mmio.c \
  ./src/libraries/matrix_reading.c \
  -o del2_r
```

**Notes:**
- `-O2`: Less aggressive optimization
- `-g`, `-Wall`, `-Wextra`: Additional debugging information

### Download Benchmark Matrices

Must be run from `Deliverable2` folder

```bash
# Download matrices from Suitsparse
./matrix_download.sh MATRIX_NAME
```

---

## Running Executables

### 1. MPI Exectuion + Sparse Matrix Generation

Run distributed sparse matrix-vector multiplication, while generating a 9% sparsity matrix of chosen size:

```bash
mpirun -np <num_ranks> ./del2_r <matrix_file> <iterations> <plot_result_file>
mpirun -np <num_ranks> ./del2_g  <iterations> <plot_result_file> <n_rows> <n_columns>
```

**Examples:**
```bash
# Run with 4 MPI working processes, 10 iterations, 256X256 matrix
mpirun -np 5 ./del2_g 10 results/to_plot/del2_g.txt 256 256

# Run with 2 MPI ranks (minimum), of which 1 working, for testing
mpirun -np 2 ./del2_g 1 results/to_plot/test.txt 256 256
```

**Output:**
- Number of Working Processes and Matrix Dimensions
- Workflow of the iterations
- Execution times (milliseconds)
- Speedup achieved

**Used for:** Weak scaling testing

**Notes:**
- Must be run with at least 2 MPI processes, since rank=0 is always the "master".
- The results in `plot_result_file` are simplified and without the workflow, to be plotted with .py scripts.

---

### 2. MPI Exectuion + Reading Matrix from File

Run distributed sparse matrix-vector multiplication, while reading matrix from file:

```bash
mpirun -np <num_ranks> ./del2_r <matrix_file> <iterations> <plot_result_file>
```

**Examples:**
```bash
# Run with 4 MPI working processes, 10 iterations
mpirun -np 5 ./del2_r matrices/11k_0p35.mtx 10 results/to_plot/del2_r.txt

# Run with 2 MPI ranks (minimum), of which 1 working, for testing
mpirun -np 2 ./del2_r matrices/11k_0p35.mtx 1 results/to_plot/test.txt
```

**Output:**
- Number of Working Processes and Matrix Dimensions
- Workflow of the iterations
- Execution times (milliseconds)
- Speedup achieved

**Used for:** Strong scaling testing

**Notes:**
- Must be run with at least 2 MPI processes, since rank=0 is always the "master".
- The results in `plot_result_file` are simplified and without the workflow, to be plotted with .py scripts.

---

### Cluster Execution (PBS)

```bash
# Submit manual matrix generation SpMV test
qsub -q short_cpuQ scripts/del2_generating.pbs

## View output
cat del2_g.out
cat del2_g.err

# Submit manual matrix reading SpMV test (example of a matrix)
qsub -q short_cpuQ -v MATRIX_FILE="matrices/11k_0p35.mtx" scripts/del2_reading.pbs

## View output
cat del2_r.out
cat del2_r.err

# Submit weak scaling test
qsub -q short_cpuQ scripts/del2_weak_scaling.pbs

## View output
cat del2_ws.out
cat del2_ws.err

# Submit strong scaling test (example of a matrix)
qsub -q short_cpuQ -v MATRIX_FILE="matrices/11k_0p35.mtx" scripts/del2_strong_scaling.pbs

## View output
cat del2_ss.out
cat del2_ss.err

# Check job status
qstat <name.username>
```

---

## Data Visualization

### Weak Scaling Visualizations

Plots the speedup in relation to the Number of Processes, separating the results based on the Matrix Size,
and the speedup in relation to the Matrix Size, separating the results based on the Number of Processes.

```bash
# Plotting of weak scaling results (results/to_plot/del2_ws.txt)
python3 scripts/plotter_ws.py
```

### Strong Scaling Visualization

Plots the speedup only in relation to the Number of Processes, since the matrix size is fixed.

```bash
# Plotting of strong scaling results (results/to_plot/del2_ss.txt)
python3 scripts/plotter_ss.py
```

**Output:** PNG files in `plots/` directory

---

## Reproducibility

**Hardware Requirements:**
- Linux cluster with MPI support
- Multi-node capability (2-128 MPI ranks)
- ≥ 32GB RAM per node

**Software Requirements:**
- MPICH 3.2.1 or compatible MPI implementation
- GCC ≥ 9.1
- Python 3 (optional): matplotlib, pandas, numpy for analysis

**Repository:** https://github.com/Suge42/PARCO-Computing-2026-242734

**Complete Workflow:**
```bash
# 1. Clone repository
git clone https://github.com/Suge42/PARCO-Computing-2026-242734.git
cd PARCO-Computing-2026-242734/Deliverable2

# 2. Download desired matrix
./scripts/matrix_download.sh desired_matrix

# 3. Cluster submission (PBS)
qsub -q short_cpuQ scripts/del2_weak_scaling.pbs
qsub -q short_cpuQ -v MATRIX_FILE="matrices/example_matrix.mtx" scripts/del2_strong_scaling.pbs

# 4. Plot the results
python3 scripts/plotter_ws.py
python3 scripts/plotter_ss.py
```

**Expected Results:**
- `results/del2_ws.out`: Weak Scaling workflow and data
- `results/to_plot/del2_ws.txt`: Weak Scaling data
- `results/del2_ss.out`: Strong Scaling workflow and data
- `results/to_plot/del2_ss.txt`: Strong Scaling data
- `results/*.err`: Eventual errors

---