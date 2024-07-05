import argparse

import numpy as np
from ase.io import read

from fcsio.reader import ForceConstantParser
from taylor import TaylorExpansionPotential


def get_original_structure(fname_in, format):
    """
    Read the original structure from a file.

    Parameters:
    - fname_in (str): The file name of the structure file.
    - format (str): The format of the structure file. If None, the format is auto-detected.

    Returns:
    - structure: The structure object read from the file.
    """
    if format is None:
        structure = read(fname_in)
    else:
        structure = read(fname_in, format=format)

    return structure


def print_energies(energies):
    """
    Print the calculated potential energies.

    Parameters:
    - energies (dict): A dictionary containing arrays of energies. The key 'total' should contain the total energy array.
    """
    print("# potential energies (Ry)")
    print("#", end='')

    keys = list(energies.keys())
    key_list = ['total'] + list(keys[:-1])

    for i, key in enumerate(key_list):
        if i == 0:
            print("{:>14s}".format(key), end='')
        else:
            print("{:>15s}".format(key), end='')
    print('')

    energies_array = np.zeros((energies['total'].shape[0], len(key_list)), dtype=float)

    for i, key in enumerate(key_list):
        energies_array[:, i] = energies[key]

    for i in range(energies['total'].shape[0]):
        for j in range(len(key_list)):
            print("{:15.7f}".format(energies_array[i, j]), end='')
        print('')


def get_argoptions():
    """
    Parse and return command line arguments.

    Returns:
    - args: The parsed command line arguments.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument('--fcs',
                        type=str, default=None,
                        help="Force constant file")

    parser.add_argument('--structure',
                        metavar='SPOSCAR', default=None,
                        help="The structure file that contains the information of original supercell")

    parser.add_argument('--format',
                        metavar='vasp', default=None,
                        help="The format of the structure specified by --structure.")

    parser.add_argument('--maxorder', type=int, default=4,
                        help="The maximum order of the Taylor expansion. (2: Harmonic, 3: Cubic, etc.)")

    parser.add_argument('--disp', type=str, default=None,
                        help="File containing the displacements in units of Bohr")

    return parser.parse_args()


def main():
    """
    Main function to execute the Taylor expansion potential calculation.
    """
    args = get_argoptions()
    structure = get_original_structure(args.structure, args.format)

    parser = ForceConstantParser(args.fcs)
    fcs_dic = parser.read(maxorder=args.maxorder)
    taylor = TaylorExpansionPotential(structure, maxorder=args.maxorder)
    taylor.set_forceconstants(fcs_dic, parser.primitive_cell)

    if args.disp is None:
        raise RuntimeError("Displacement file is not specified.")

    displacements = np.loadtxt(args.disp)
    nlen = len(displacements) // (structure.get_global_number_of_atoms())
    if nlen * structure.get_global_number_of_atoms() != len(displacements):
        raise RuntimeError("The number of displacements is inconsistent.")

    displacements = displacements.reshape(-1, structure.get_global_number_of_atoms(), 3)
    energies, forces = taylor.compute(displacements)
    print_energies(energies)


if __name__ == '__main__':
    main()
