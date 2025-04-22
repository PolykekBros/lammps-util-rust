import subprocess
import sys
from pathlib import Path
import math

def main():
    dir = Path(sys.argv[1])
    cnt = 0
    sum = [0, 0, 0]
    sum2 = [0, 0, 0]
    for i in range(0,100):
        run_dir = dir / f"run_{i+1}"
        initial = run_dir / "dump.initial"
        final = run_dir / "dump.final_no_cluster"
        # print(initial, final)
        output = subprocess.run(
            [
                "../target/release/component-shift",
                initial,
                final,
            ],
            capture_output=True, check=True, text=True
        ).stdout
        lines = str(output).split('\n')
        cnt += float(lines[0])
        sum_tmp = [float(n) for n in lines[1].split()]
        sum2_tmp = [float(n) for n in lines[2].split()]
        for i in range(0,3):
            sum[i] += sum_tmp[i]
            sum2[i] += sum2_tmp[i]
        # print(cnt, sum_tmp, sum2_tmp)
    print(cnt, sum, sum2)
    print(sum[0] / cnt, sum[1] / cnt, sum[2] / cnt)
    print(math.sqrt(sum2[0] / cnt), math.sqrt(sum2[1] / cnt), math.sqrt(sum2[2] / cnt))

    
if __name__ == "__main__":
    main()
