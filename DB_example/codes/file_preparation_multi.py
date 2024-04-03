import mdtraj as md
import numpy as np
import sys
import csv


#Recognise each argument
for i in range(0, len(sys.argv)):
    if sys. argv[i] == "--protein_ligand":
        protein_ligand = sys.argv[i + 1]
    elif sys.argv[i] == "--replicas":
        replicas = int(sys.argv[i + 1])
    if len(sys.argv) < 5:
        print('Error: python file_preparation.py --protein-ligand input.csv --replicas n_replicas')
        print('input.csv: complex_name, chainA, chainB, SMILES')
        exit()


################################################
#                                              #
#                READ CSV FILE                 #
#                                              #
################################################
def read_csv(protein_ligand):
    # Open the CSV file
    with open(protein_ligand, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Skip the header row
        data = []
        for row in reader:
            complex_name, chainA, chainB, SMILES = row
            data.append((complex_name, chainA, chainB, SMILES))
        return data


################################################
#                                              #
#               MERGE PDB FILES                # Asegurarse que las estructuras estan separadas
#                                              #
################################################
def distance(xyz1, xyz2):

    '''
    This function ensures that the structures are separated by at least 10 angstroms.

    Inputs:
    xyz1 - The coordinates of the first structure.
    xyz2 - The coordinates of the second structure.

    Output:
    xyz2 - The coordinates of the second structure, shifted by 10 angstroms if necessary. 
      
    '''

    if np.linalg.norm(xyz1.mean(axis = 1) - xyz2.mean(axis = 1)) < 10:
        mean1 = xyz1.mean(axis = 1)
        mean2 = xyz2.mean(axis = 1)
        vector = mean2 - mean1
        distance = np.linalg.norm(vector)
        unit_vector = vector / distance
        separation_vector = unit_vector * 10
        xyz2 += separation_vector


################################################
#                                              #
#               CREATE .csv FILE               # Crear el archivo CSV que será en input de Diffdock
#                                              #
################################################
def create_csv(data):

    '''
    Creates a .csv file with the protein path and SMILES.

    Imputs:
    protein_path - A string representing the path to the protein. (e.g. '1a07.pdb') si este archivo esta en pdb_files
    ligand_description - A string representing the SMILES of the ligand. (e.g. 'CNCC1=CC=C(C=C1)C2=C3CCNC(=O)C4=C3C(=CC(=C4)F)N2')

    Output:
    protein_ligand.csv - A .csv file with the columns: complex_name, protein_path, ligand_description, protein_sequence.   
    
    '''

    with open('protein_ligand.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["complex_name", "protein_path", "ligand_description", "protein_sequence"])
        for i in range(replicas):
            for row in data:
                complex_name, chainA, chainB, SMILES = row
                writer.writerow([complex_name, complex_name + '.pdb', SMILES, ""])

# Main code
if __name__ == '__main__':
    data = read_csv(protein_ligand)

    for complex_name, chainA, chainB, SMILES in data:
        # Now you can use these variables in your script

        #Load the structures
        p1 = md.load(chainA)
        p2 = md.load(chainB)

        top1, xyz1 = p1.top, p1.xyz
        top2, xyz2 = p2.top, p2.xyz
        top1_2 = top1.join(top2)

        distance(xyz1, xyz2)

        #Merge the structures and save the new pdb file
        xyz1_2 = np.concatenate((xyz1, xyz2), axis=1)
        p1_2 = md.Trajectory(xyz=xyz1_2, topology=top1_2) 
        p1_2.save(complex_name + '.pdb')
        print('New pdb file created: ' + complex_name + '.pdb')

        file = complex_name + '.pdb' # es lo mismo que se pone en --out del merge.py pero con .pdb
        #SMILES = 'CNCC1=CC=C(C=C1)C2=C3CCNC(=O)C4=C3C(=CC(=C4)F)N2'
        create_csv(data)

    print('New .csv file created: protein_ligand.csv')
   

