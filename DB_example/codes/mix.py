import numpy as np
import pandas as pd
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument(
    "--data_a",
    type=str,
    help="Path to the file containing the proteins and ligands with their affinities of proteins A (more affine). Must be csv",
    required=True,
)
parser.add_argument(
    "--data_b",
    type=str,
    help="Path to the file containing the proteins and ligands with their affinities of proteins B (less affine). Must be csv",
    required=True,
)
parser.add_argument(
    "--th_amb",
    type=float,
    help="Threshold from which the data will be considered ambiguous. Default is 0.5",
    default=0.5,
    required=False,
)
args = parser.parse_args()

th = args.th_amb # Ka/Kb that will be considered as "equal" and therefore non differentiable

# ensure that the Training_data dir exists

results = os.listdir("DD")

dfa = pd.read_csv(args.data_a)
dfb = pd.read_csv(args.data_b)

# select the values that are ambiguous
mask = dfa[dfa["Kd (nM)"] / dfb["Kd (nM)"] >= th] 
dfa = dfa.drop(mask.index)
dfb = dfb.drop(mask.index)

id_c = []
for i in range(len(dfa)):
    id_c.append(
        f"{dfa['Prot ID'].values[i]}_{dfb['Prot ID'].values[i]}_{int(dfa['PubChem CID'].values[i])}"
    )


for r in results:
    df = pd.read_csv("DD/" + r)
    mask = df["ID"].isin(id_c) # Separate the ambiguous values.
    dfc = df[mask].copy()
    total = np.unique(df[~mask]["ID"]) # Now the total IDs doesn't contain the ambiguous values
    select_1 = np.random.choice(total, size=(len(total) // 2)) # Half the values will change order
    select_0 = np.setdiff1d(total, select_1) # The other half will remain unchanged

    df_1 = df[df["ID"].isin(select_1)].copy()
    df_0 = df[df["ID"].isin(select_0)].copy()

    ind_1 = df_1.index
    ind_0 = df_0.index

    # Swap prot A and prot B for the selected rows.
    id = df_1["ID"].values
    id = [f"{i.split('_')[1]}_{i.split('_')[0]}_{i.split('_')[2]}" for i in id]
    df_1["ID"] = id

    # Add a colum with the experimental prediction:
        # A>B is represented as 'A'
        # B>A is represented as 'B'
        # Ambiguity is represented as 'C'  
    df_1["Experimental Pred"] = "B"
    df_0["Experimental Pred"] = "A"
    dfc["Experimental Pred"] = "C"

    # Change 0 for 1 on the selected chains.
    df_1.loc[df_1["Chain(A=0)(B=1)"] == 0, "Chain(A=0)(B=1)"] = -1
    df_1.loc[df_1["Chain(A=0)(B=1)"] == 1, "Chain(A=0)(B=1)"] = 0
    df_1.loc[df_1["Chain(A=0)(B=1)"] == -1, "Chain(A=0)(B=1)"] = 1

    # Put the three kind of values together.
    df_1 = pd.concat([df_1 if not df_1.empty else None, df_0], ignore_index=True)
    df_1 = pd.concat([df_1 if not df_1.empty else None, dfc], ignore_index=True)

    df_1.to_csv("model/training_data/" + r.split('.')[0] + '_train.csv', index=False)

print('')
print('Training data successfully generated')