import pandas as pd
import numpy as np
import os
from itertools import combinations
import mdtraj as md
from argparse import ArgumentParser

def args_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "--SMILES",
        type=str,
        help="Canonical SMILES of the ligand.",
        required=True,
    )
    parser.add_argument(
        "--run_name",
        type=str,
        help="Name of the job.",
        required=True,
    )
    args = parser.parse_args()
    return args

def new_pdb(output_pdb, prot_a, prot_b, separate_chains=[10,0,0]):
    """
    This function concatenates multiple mdtraj objects into a single PDB file spacing them by a factor.

    Parameters
    ----------
    output_pdb : str. The path to the output PDB file.
    chains : list. A list of mdtraj objects.
    separate_chains : (list, optional) The distance (in Angstroms) to separate each chain, by default [10, 0, 0]

    Returns
    -------
    None. Generates a pdb file with the merged proteins properly spaced.

    """

    prot_a = md.load(prot_a)
    prot_b = md.load(prot_b)

    top_a = prot_a.top # extract protein topology
    top_b = prot_b.top 

    xyz_a = prot_a.xyz # extract the protein coordinates
    xyz_b = prot_b.xyz
 
    xyz_b += separate_chains # apply separation desired for every chain

    merged_top = top_a.join(top_b) # merge into a single pdb
    merged_xyz = np.concatenate((xyz_a, xyz_b),axis=1)

    merged_pdb = md.Trajectory(xyz=merged_xyz, topology=merged_top) 

    merged_pdb.save(f"{output_pdb}")
    print(f'File {output_pdb} generated successfully')

def generate_input(structure_dir, smiles, run_name):

    # Get all the structures
    structures = os.listdir(structure_dir)
    
    try:
        os.mkdir(f'proteins/paired_structures/{run_name}')
    except:
        pass
    
    n_structures = len(structures)

    # Generate all possible pairs
    ind = np.arange(n_structures)
    pair_ind = list(combinations(ind, 2))

    print(f'--> Using {n_structures} structures to generate {len(pair_ind)} different possible pairs.')

    # Generate the input file
    input = pd.DataFrame(columns=['complex_name','ligand_description','protein_path','protein_sequence'])
    ids = []
    s = []
    p = []


    for pair in pair_ind:
        prot_a = structures[pair[0]]
        prot_b = structures[pair[1]]
        id = f"{prot_a.split('.')[0]}_{prot_b.split('.')[0]}"

        new_pdb(f'proteins/paired_structures/{run_name}/{id}.pdb', f'{structure_dir}/{prot_a}', f'{structure_dir}/{prot_b}')

        ids.append(id)
        s.append(smiles)
        p.append(os.path.abspath(f'proteins/paired_structures/{run_name}/{id}.pdb'))

    input['complex_name'] = ids
    input['ligand_description'] = s
    input['protein_path'] = p

    input.to_csv(f'inputs/input_{run_name}.csv', index=False)

    



if __name__ == '__main__':
    args = args_parser()
    print('')
    print(f'Generating input for {args.run_name}')
    generate_input(f'proteins/protein_structures/{args.run_name}',args.SMILES,args.run_name)