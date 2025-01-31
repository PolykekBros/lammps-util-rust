import math
import sys
import subprocess

target = sys.argv[1]
sim_n = 100
cutoff_start = 1.6
cutoff_end = 2.0
cutoff_step = 0.05

cutoffs_len = int(math.ceil((cutoff_end - cutoff_start) / cutoff_step)) + 1
cutoffs = [cutoff_start + i * cutoff_step for i in range(cutoffs_len)]
print(f"cutoffs: {cutoffs}")


def process(s):
    return [float(t) for t in s.strip().split()]


def run(c):
    avg = [0, 0, 0, 0, 0]
    for i in range(sim_n):
        dir = f"{target}/run_{i + 1}"
        di = f"{dir}/dump.initial"
        df = f"{dir}/dump.final_no_cluster"
        out = subprocess.run(
            ["crater-analysis", di, df, dir, "-c", str(c)],
            check=True,
            capture_output=True,
            text=True,
        )
        res = process(out.stdout)
        for i, r in enumerate(res):
            avg[i] += r
    for i in range(len(avg)):
        avg[i] /= sim_n
    return avg


for c in cutoffs:
    avg = run(c)
    for n in [c] + avg:
        print(f"{n:8.2f}", end="")
    print()
