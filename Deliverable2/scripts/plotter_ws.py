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
df = parse_benchmark_file("results/to_plot/del2_ws.txt")

# Add a convenience column
df["matrix_size"] = df["rows"] * df["cols"]

#print(df)



# Plotting speedup vs number of processes
plt.figure(figsize=(8, 5))

## One curve per matrix size
for size, group in df.groupby("matrix_size"):
    group = group.sort_values("processes")
    plt.plot(
        group["processes"],
        group["speedup"],
        marker="o",
        label=f"{int(np.sqrt(size))}×{int(np.sqrt(size))}"
    )

## Ideal linear speedup
max_p = df["processes"].max()
plt.plot(
    range(1, max_p + 1),
    range(1, max_p + 1),
    "k--",
    label="Ideal speedup"
)

plt.xlabel("Number of processes")
plt.ylabel("Speedup")
plt.title("Speedup vs Number of Processes")
plt.legend(title="Matrix size")
plt.tight_layout()

plt.savefig("plots/del2_ws_speedup_vs_processes.png", dpi=300)
plt.close()




# Plotting speedup vs matrix size
plt.figure(figsize=(8, 5))

## One curve per number of processes
for p, group in df.groupby("processes"):
    group = group.sort_values("matrix_size")
    plt.plot(
        group["matrix_size"],
        group["speedup"],
        marker="o",
        label=f"{p} processes"
    )

plt.xscale("log")
plt.xlabel("Matrix size (rows × cols)")
plt.ylabel("Speedup")
plt.title("Speedup vs Matrix Size")
plt.legend(title="Parallelism")
plt.tight_layout()

plt.savefig("plots/del2_ws_speedup_vs_matrix_size.png", dpi=300)
plt.close()
