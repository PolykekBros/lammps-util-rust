import sys

target = sys.argv[1]
sim_n = 100


def process(s):
    return [float(t) for t in s.strip().split()[1:]]


avg = [0, 0, 0, 0, 0]


with open(target) as f:
    for line in f.readlines():
        res = process(line)
        for i, r in enumerate(res):
            avg[i] += r

for i in range(len(avg)):
    avg[i] /= sim_n

for n in avg:
    print(f"{n:8.2f}", end="")
print()
