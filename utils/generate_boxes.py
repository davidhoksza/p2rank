"""
Generates pocket boxes ready to be used for docking.
"""
import sys
import argparse
import logging
import glob
from Bio.PDB import *
import numpy as np
import json

from typing import List, Dict


class Pocket:
    def __init__(self, id: str, center: List[float], residues: List[str], atoms: List[str]):
        self.id = id
        self.center: List[float] = center
        self.residues: List[str] = residues
        self.atoms: List[str] = atoms
        self.bb: List = None #box size

def error_exit(message):
    logging.error(message)
    sys.exit()


def parse_structure(path: str):
    parser = PDBParser()
    return parser.get_structure("str", path)


def parse_pocket(line: str) -> Pocket:
    s_line = line.split(',')
    return Pocket(
        id=s_line[0],
        center=[float(s_line[5]), float(s_line[6]), float(s_line[7])],
        residues=[r.strip() for r in s_line[8].split(' ')],
        atoms=[a.strip() for a in s_line[9].split(' ')],)


def parse_pockets(path: str) -> List[Pocket]:
    pockets: List[Pocket] = []
    with open(path, 'r') as f:
        for line in list(f.readlines())[1:]:
            pockets.append(parse_pocket(line))
    return pockets


def get_all_atoms(structure: 'Bio.PDB.Structure') -> Dict[str, 'Bio.PDB.Atom']:
    atoms: Dict[str, 'Bio.PDB.Atom'] = {}
    for chain in structure[0]:
        for residue in chain:
            for atom in residue:
                atoms[atom.get_serial_number()] = atom
    return atoms


def get_atom_coords(atom_id: str, atoms: Dict[str, 'Bio.PDB.Atom']) -> List[float]:
    return atoms[atom_id].get_coord()


def attach_bb(pocket: Pocket, atoms: Dict[str, 'Bio.PDB.Atom'], padding: float):
    max_diff = np.array([0,0,0])
    center = np.array(pocket.center)
    for atom_id in pocket.atoms:
        coords = np.array(atoms[int(atom_id)].get_coord())
        diff = np.abs(coords-center)
        max_diff = np.maximum(diff, max_diff)

    pocket.bb = list(max_diff + [padding, padding, padding])


def attach_bbs(pockets: List[Pocket], structure: 'Bio.PDB.Structure', padding: float):
    atoms = get_all_atoms(structure=structure)
    for pocket in pockets:
        attach_bb(pocket, atoms, padding)


def get_pockets_with_boxes(dir: str, padding: float) -> List[Pocket]:

    pockets_path: str = list(glob.iglob('{}/*_predictions.csv'.format(dir)))[0]
    pockets: List[Pocket] = parse_pockets(pockets_path)

    pdb_path:str = list(glob.iglob('{}/visualizations/data/*.pdb'.format(dir)))[0]
    structure = parse_structure(pdb_path)

    attach_bbs(pockets, structure, padding)

    return pockets


def main():

    pockets = get_pockets_with_boxes(dir=args.input, padding=float(args.padding))
    with (sys.stdout if args.output is None else open(args.output, "w")) as fw:
        fw.write(json.dumps([{"id": p.id,
                              "center": p.center,
                              "box": p.bb
                              } for p in pockets], indent=4))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input",
                        required=True,
                        metavar='DIR',
                        help="The path with the P2Rank output")
    parser.add_argument("-p", "--padding",
                        metavar='FLOAT',
                        default=0,
                        help="Number of Angstromes to be added to each dimension of the minimum bounding box")
    parser.add_argument("-o", "--output",
                        metavar='FILE',
                        help="Output file name. If non entered, the standard output will be used.")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')

    main()