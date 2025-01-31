import math
import sys
import asyncio
from asyncio import subprocess

target = sys.argv[1]
sim_n = 100
cutoff_start = 1.745
cutoff_end = 1.755
cutoff_step = 0.005

n_proc = 10

cutoffs_len = int(math.ceil((cutoff_end - cutoff_start) / cutoff_step)) + 1
cutoffs = [cutoff_start + i * cutoff_step for i in range(cutoffs_len)]
print(f"cutoffs: {cutoffs}")


async def exec(*args):
    proc = await subprocess.create_subprocess_exec(
        *args,
        stdout=subprocess.PIPE,
    )
    code = await proc.wait()
    if code != 0:
        sys.exit(1)
    stdout, _ = await proc.communicate()
    return stdout.decode()


def process(s):
    return [float(t) for t in s.strip().split()]


async def run_one(i, c):
    dir = f"{target}/run_{i + 1}"
    di = f"{dir}/dump.initial"
    df = f"{dir}/dump.final_no_cluster"
    out = await exec("crater-analysis", di, df, dir, "-c", str(c))
    return process(out)


async def run_n(i_start, n, c):
    tasks = [run_one(i_start + i, c) for i in range(n)]
    return await asyncio.gather(*tasks)


def run(c):
    avg = [0, 0, 0, 0, 0]
    for i in range(sim_n // n_proc):
        res = asyncio.run(run_n(i * n_proc, n_proc, c))
        for row in res:
            for i, r in enumerate(row):
                avg[i] += r
    for i in range(len(avg)):
        avg[i] /= sim_n
    return avg


for c in cutoffs:
    avg = run(c)
    for n in [c] + avg:
        print(f"{n:9.3f}", end="")
    print()
