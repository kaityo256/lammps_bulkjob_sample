import subprocess


def main():
    temp_array = [0.9, 0.92, 0.94, 0.98, 1.00, 1.02, 1.04, 1.06]
    density = 0.9
    for t in temp_array:
        cmd = ["python3", "generate_config.py", "-d", str(density), "-t", str(t)]
        print(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
