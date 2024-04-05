import numpy as np
import mdtraj as md
import os
import re
import pandas as pd


def get_confidence(dir):
    """
    This function takes a directory as input and returns a list of the confidence values of the models in the directory, sorted in descending order.

    Parameters:
    dir (str): Path to the directory containing the result files of the DifDock simulation (converted from sdf to pdb).

    Returns:
    list: A list of the confidence values of the models, sorted in descending order.

    """
    conf_list = []
    files = os.listdir(dir)
    for file in files:
        # print(file)
        try:
            conf = re.search(r"confidence(.*?)\.pdb", file)
            conf = conf.group(1)
            conf_list.append(float(conf))
        except AttributeError:
            continue

    sorted_list = sorted(conf_list, reverse=True)
    sorted_list = [f"{i:.2f}" for i in sorted_list]
    assert (
        len(sorted_list) != 0
    ), f"None of the files in {dir} can be analyzed... maybe they are in sdf format instead of pdb?"
    return sorted_list


def lig_gc(lig_xyz):
    """
    This function takes a numpy array of coordinates as input and returns the geometric center of the ligand as a numpy array.

    Parameters:
    lig_xyz (np.ndarray): A numpy array of shape (n_frames, n_atoms, 3) containing the coordinates of the ligand in Cartesian coordinates.

    Returns:
    np.ndarray: A numpy array of shape (3,) containing the geometric center of the ligand.

    """

    a = np.array(lig_xyz[0])
    mean = np.mean(a, axis=0)
    gc = [mean[0], mean[1], mean[2]]  # x, y, z
    return gc


def get_std_prot(dir, th=0.2):

    conf = get_confidence(dir)
    std_sim_old = 0
    gc_list = []
    
    for n, c in enumerate(conf):
        
        # avoid problems with confidence = 0.0 or -0.0
        if abs(float(c)) == 0:
            try:
                file = f"rank{n+1}_confidence0.00.pdb"
                lig = md.load(dir + file)
            except OSError:
                file = f"rank{n+1}_confidence-0.00.pdb"
                lig = md.load(dir + file)

        else:
            file = f"rank{n+1}_confidence{c}.pdb"
            lig = md.load(dir + file)

        lig_xyz = lig.xyz
        if np.isnan(lig_xyz).any():
            # print(f'Coordinates of {file} may be infinite, skipping file in std calculation')
            pass
        else:
            gc = lig_gc(lig_xyz)
            gc_list.append(gc)
            std_sim_new = np.linalg.norm(np.std(np.array(gc_list), axis=0))
            if std_sim_new == 0:
                continue
            if abs(std_sim_old - std_sim_new) < 2:  # std difference threshold
                std_sim_old = std_sim_new
            else:
                # print(std_sim_old, n+1)
                if std_sim_old > th:
                    fail = 1
                    # print('DiffDock may have *FAILED* to predict the docking')
                else:
                    # print('DiffDock seems to have *CORRECTLY* predicted the docking')
                    fail = 0
                return std_sim_old, fail

    # print(std_sim_old,n+1)
    if std_sim_old > th:
        fail = 1
        # print('DiffDock may have *FAILED* to predict the docking site')
    else:
        # print('DiffDock seems to have *CORRECTLY* predicted the docking site')
        fail = 0
    return std_sim_old, fail


if __name__ == "__main__":
    fails = {'Fail':[],'Protein':[],'Pubchem CID':[]}
    x_std=[]
    # n_samples = np.arange(0,n+1)

    for dir_path in os.listdir('/home/ramon/juan/DD/DD_solo/'):
         if os.path.isdir(f'/home/ramon/juan/DD/DD_solo/{dir_path}'):
            s, fail = get_std_prot(f'/home/ramon/juan/DD/DD_solo/{dir_path}/')
            name = dir_path.split('_')
            fails['Fail'].append(fail)
            fails['Protein'].append(name[0])
            fails['PubChem CID'].append(name[1])
    df = pd.DataFrame(fails)
    df.to_csv('fails.csv',index=False)