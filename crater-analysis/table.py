import asyncio
from asyncio import subprocess
import sys

target = sys.argv[1]
sim_n = 100
cutoff = 1.75

n_proc = 10


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


async def run_one(i, c):
    i = i + 1
    dir = f"{target}/run_{i}"
    di = f"{dir}/dump.initial"
    df = f"{dir}/dump.final_no_cluster"
    return (i, await exec("crater-analysis", di, df, dir, "-c", str(c)))


async def run_n(i_start, n, c):
    tasks = [run_one(i_start + i, c) for i in range(n)]
    return await asyncio.gather(*tasks)


lines = []
for i in range(sim_n // n_proc):
    lines += asyncio.run(run_n(i * n_proc, n_proc, cutoff))

lines = sorted(lines, key=lambda p: p[0])
lines = [f"{p[0]} {p[1]}" for p in lines]
print("".join(lines), end="")
