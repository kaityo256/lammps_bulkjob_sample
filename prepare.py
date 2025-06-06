import subprocess


def save_jobfile(temp_array):
    str = ""
    for temperature in temp_array:
        prefix = f"T{int(temperature*100):03d}"
        str += f"srun --exclusive --mem-per-cpu=1840 -n 128 -c 1 lammps < {prefix}/{prefix}.input &\n"
        str += f"sleep 5\n"
    content = f"""\
#!/bin/sh

#SBATCH -p i8cpu
#SBATCH -N 8

source /home/issp/materiapps/intel/lammps/lammpsvars.sh
{str}
wait
"""
    filename = "data/job8.sh"
    with open("data/job8.sh", "w") as f:
        f.write(content)
    print(f"Generated {filename}")


def main():
    temp_array = [0.9, 0.92, 0.94, 0.96, 0.98, 1.00, 1.02, 1.04]
    density = 0.9
    save_jobfile(temp_array)
    for t in temp_array:
        cmd = ["python3", "generate_config.py", "-d", str(density), "-t", str(t)]
        print(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
