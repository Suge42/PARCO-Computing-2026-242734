import os
import matplotlib.pyplot as plt
import re
import pandas as pd
import numpy as np

# Output directory for
os.makedirs("plots", exist_ok=True)
plt.style.use("seaborn-v0_8-darkgrid")


def parse_benchmark_file(filename):
    data = []
    current = {}

    matrix_pattern = re.compile(
        r"Row:\s*(\d+)\s*-\s*Columns:\s*(\d+)\s*-\s*Processes:\s*(\d+)"
    )
    value_pattern = re.compile(r"(\w+):\s*([0-9.]+)")

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            if line.startswith("#"):
                match = matrix_pattern.search(line)
                if match:
                    if "speedup" in current:
                        data.append(current)
                        current = {}

                    current["rows"] = int(match.group(1))
                    current["cols"] = int(match.group(2))
                    current["processes"] = int(match.group(3))
            else:
                match = value_pattern.match(line)
                if match:
                    current[match.group(1)] = float(match.group(2))

        if current:
            data.append(current)

    return pd.DataFrame(data)




# Read and parse the benchmark file
df = parse_benchmark_file("results/to_plot/del2_ss.txt")

# Add a convenience column
df["matrix_size"] = df["rows"] * df["cols"]

# Sanity check: ensure single matrix size
matrix_sizes = df["matrix_size"].unique()
if len(matrix_sizes) != 1:
    raise ValueError(
        f"Expected exactly one matrix size, found {len(matrix_sizes)}: {matrix_sizes}"
    )

matrix_size = matrix_sizes[0]
n = int(np.sqrt(matrix_size))

# Sort by number of processes
df = df.sort_values("processes")

# Plotting speedup vs number of processes
plt.figure(figsize=(8, 5))

plt.plot(
    df["processes"],
    df["speedup"],
    marker="o",
    linewidth=2,
    label="Measured speedup"
)

# Ideal linear speedup
max_p = df["processes"].max()
plt.plot(
    range(1, max_p + 1),
    range(1, max_p + 1),
    "k--",
    label="Ideal speedup"
)

plt.xlabel("Number of processes")
plt.ylabel("Speedup")
plt.title(f"Speedup vs Number of Processes ({n}×{n} matrix)")
plt.legend()
plt.tight_layout()

plt.savefig("plots/del2_ss_speedup_vs_processes.png", dpi=300)
plt.close()
