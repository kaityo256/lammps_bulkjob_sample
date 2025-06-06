import argparse
import random
import sys
import os


class Atom:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.type = 1
        V0 = 1.0
        self.vx = V0 * random.uniform(-1, 1)
        self.vy = V0 * random.uniform(-1, 1)
        self.vz = V0 * random.uniform(-1, 1)


def make_config(m, s):
    s = 1.55
    h = 0.5 * s
    atoms = []
    for ix in range(0, m):
        for iy in range(0, m):
            for iz in range(0, m):
                x = ix * s
                y = iy * s
                z = iz * s
                atoms.append(Atom(x, y, z))
                atoms.append(Atom(x, y + h, z + h))
                atoms.append(Atom(x + h, y, z + h))
                atoms.append(Atom(x + h, y + h, z))
    return atoms


def save_atoms_file(dirname, atoms, L):
    filename = f"{dirname}/input.atoms"
    with open(filename, "w") as f:
        f.write("Position Data\n\n")
        f.write(f"{len(atoms)} atoms\n")
        f.write("1 atom types\n\n")
        f.write(f"0 {L} xlo xhi\n")
        f.write(f"0 {L} ylo yhi\n")
        f.write(f"0 {L} zlo zhi\n")
        f.write("\n")
        f.write("Atoms\n\n")
        for i, a in enumerate(atoms):
            f.write(f"{i+1} {a.type} {a.x} {a.y} {a.z}\n")
        f.write("\n")
        f.write("Velocities\n\n")
        for i, a in enumerate(atoms):
            f.write(f"{i+1} {a.vx} {a.vy} {a.vz}\n")
    print(f"Generated {filename}")


def save_input_file(dirname, temperature):
    prefix = f"T{int(temperature*100):03d}"
    filename = f"{dirname}/{prefix}.input"
    content = f"""\
units lj
atom_style atomic
boundary p p p
timestep 0.001

read_data input.atoms

mass 1 1.0

pair_style lj/cut 3.0
pair_coeff 1 1 1.0 1.0 3.0

fix 1 all nvt temp {temperature} {temperature} 1
log {prefix}.log
thermo 100
run 1000
"""
    with open(filename, "w") as f:
        f.write(content)
    print(f"Generated {filename}")


def save_files(density, temperature):
    S = (4.0 / density) ** (1.0 / 3.0)  # Lattice Constant
    M = 10  # Lattice Count
    N = M**3 * 4
    L = M * S  # Lattice Size
    atoms = make_config(M, S)
    density = N / (L**3)
    prefix = f"T{int(temperature*100):03d}"
    dirname = f"data/{prefix}"
    os.makedirs(dirname, exist_ok=True)
    save_atoms_file(dirname, atoms, L)
    save_input_file(dirname, temperature)


def main():
    parser = argparse.ArgumentParser(
        description="Specify the density and temperature.",
        usage="%(prog)s -d DENSITY -t TEMPERATURE",
    )
    parser.add_argument(
        "-d", "--density", type=float, required=True, help="Density value (float)."
    )
    parser.add_argument(
        "-t",
        "--temperature",
        type=float,
        required=True,
        help="Temperature value (float).",
    )
    # Show usage if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_usage()
        sys.exit(1)

    args = parser.parse_args()
    save_files(args.density, args.temperature)


if __name__ == "__main__":
    random.seed(0)
    main()
