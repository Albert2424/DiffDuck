import numpy as np
import pandas as pd
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument(
    "--results_file",
    type=str,
    help="Path to the file containing the DD results. Must be csv",
    required=True,
)
parser.add_argument(
    "--output_file",
    type=str,
    help="File that will be saved in output directory and will contain the predictions of the DD run.",
    required=True,
)
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
args = parser.parse_args()


def prediction(df, out_file):
    """
    This function takes a pandas dataframe as input and predicts the chain in which the ligand has attached.
    The prediction is based on the confidence scores of each prediction chain.

    Args:
        df (pd.DataFrame): A pandas dataframe containing the prediction results and customer information.
                            The columns of the dataframe must include:
                                'ID', 'Confidence', 'Chain(A=0)(B=1)'

    Returns:
        pd.DataFrame: A pandas dataframe containing the chain ID, prediction, and Reliability of the prediction.
    """


    pred_df = pd.DataFrame(columns=["ID", "Prediction", "Reliability"])

    for i in np.unique(df["ID"]):

        slice = df[df["ID"] == i].copy()
        slice = slice[slice["Confidence"] != -1000]  # exclude failed attempts
        slice = slice.reset_index(drop=True)
        # slice = slice.iloc[:50]
        samples = len(slice)

        slice.loc[slice["Chain(A=0)(B=1)"] == 0, "Chain(A=0)(B=1)"] = -1
        slice.loc[slice['Chain(A=0)(B=1)'] == 1, 'Chain(A=0)(B=1)'] = 0
        slice.loc[slice['Chain(A=0)(B=1)'] == -1, 'Chain(A=0)(B=1)'] = 1
        # slice.loc[slice['Chain(A=0)(B=1)'] == 0, 'Chain(A=0)(B=1)'] = -1

        ci = np.array(slice["Confidence"])
        c_max = ci[0]
        c_min = ci[len(slice) - 1]

        numerator = ci - c_min
        denominator = c_max - c_min

        norm_conf = np.array(numerator) / denominator

        slice["Normal conf"] = norm_conf
        d_0 = slice.loc[slice["Chain(A=0)(B=1)"] == 0, "Normal conf"] 
        c_0 = np.sum(np.array(d_0))/np.sum(norm_conf)
   
        d_1 = slice.loc[slice["Chain(A=0)(B=1)"] == 1, "Normal conf"]
        c_1 = np.sum(np.array(d_1))/np.sum(norm_conf)

        # ni/N
        # w_0 = len(d_0)/samples
        # w_1 = len(d_1)/samples

        # sum(ind/samples) / w_max 
        # w_max = np.sum(1/np.arange(1,samples))
        # w_0 = np.sum(1/(np.array(d_0.index) + 1))/w_max
        # w_1 = np.sum(1/(np.array(d_1.index) + 1))/w_max

        # simply 1
        w_0 = 1
        w_1 = 1

        # change weight every 10 samples
        # w_max = np.sum([(i+1) for i in range(int(samples/10))])
        # ind_0 = np.array(d_0.index) + 1
        # ind_1 = np.array(d_1.index) + 1

        # i_0 = []
        # i_1 = []
        # factor = 10
        # for j in (ind_0):
        #     if j < factor:
        #         i_0.append(1/factor)
        #     else:
        #         factor += 10
        #         i_0.append(1/factor)

        # factor = 10
        # for j in ind_1:
        #     if j < factor:
        #         i_1.append(1/factor)
        #     else:
        #         factor += 10
        #         i_1.append(1/factor) 

        # w_0 = np.sum((np.array(i_0)))/w_max
        # w_1 = np.sum((np.array(i_1)))/w_max
        
        pred_0 = c_0*w_0
        pred_1 = c_1*w_1

        # print(pred_0,pred_1) 
        # print(w_0,w_1, c_0, c_1)   

        if pred_1 > pred_0:
            pred_df = pd.concat(
                [
                    pred_df if not pred_df.empty else None,
                    pd.DataFrame(
                        [
                            {
                                "ID": i,
                                "Prediction": 1,
                                "Reliability": round(abs(pred_1) * 100, 2),
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )
            # print(f'The prediction is chain 1 with certainty of {pred_1*100:.2f}%')
        else:
            pred_df = pd.concat(
                [
                    pred_df if not pred_df.empty else None,
                    pd.DataFrame(
                        [
                            {
                                "ID": i,
                                "Prediction": 0,
                                "Reliability": round(abs(pred_0) * 100, 2),
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )
            # print(f'The prediction is chain 0 with certainty of {(1-pred_1)*100:.2f}%')
        # print(pred_1)

        pred_df.to_csv("output/" + out_file)
    return pred_df


def plot_rel_rat(dfa, dfb, pred_df, out):
    import matplotlib.pyplot as plt

    # if DD has failed there will be different values at the dfa and dfb so we must avoid them
    for i in range(len(dfa["Prot ID"])):
        id = f"{dfa['Prot ID'].loc[i]}_{dfb['Prot ID'].loc[i]}_{int(dfa['PubChem CID'].loc[i])}"
        # print(id)
        if id in pred_df["ID"].values:
            pass
        else:
            dfa = dfa.drop(i)
            dfb = dfb.drop(i)

    # pred_df = pred_df[pred_df['Prediction'] == 1]
    x = np.array(pred_df["Reliability"])
    y = np.log10(np.array(dfb["Kd (nM)"]) / np.array(dfa["Kd (nM)"]))

    col = np.array(pred_df["Prediction"], dtype=float)
    col[col == 0] = 0.20
    correct = len(col[col == 0.2])
    failed = len(col) - correct

    if len(y) != len(col):
        print(
            f'{out.split(".")[0]}_rel.pdf will not be generated since the prediction file and the data_A and data_B are not corresponding: {dfa.shape}, {len(col)}'
        )
        return

    plt.figure(figsize=(20, 15))
    a = plt.scatter(x, y, c=col, cmap="PiYG_r", s=35, clim=(0, 1))
    plt.xlabel("Reliability", fontsize=35)
    plt.ylabel(rf"$log(K_d^B/K_d^A)$", fontsize=35)
    plt.legend(
        fontsize=35,
        markerscale=4.0,
        handles=a.legend_elements()[0],
        labels=[f"Guessed: {correct}", f"Failed: {failed}"],
        ncol=2,
        bbox_to_anchor=(0.5, 1.05),
        fancybox=True,
        shadow=True,
        loc="upper left",
    )

    plt.xticks(fontsize=35)
    plt.yticks(fontsize=35)
    plt.tight_layout()
    plt.savefig(f"output/figures/{out.split('.')[0]}_rel.pdf")
    # plt.show()


if __name__ == "__main__":
    print(f"Ranking {args.results_file.split('/')[-1]} values...")
    df = pd.read_csv(args.results_file)
    dfa = pd.read_csv(args.data_a)
    dfb = pd.read_csv(args.data_b)

    # dfa = dfa.loc[[i for i in range(51)]]
    # dfb = dfb.loc[[i for i in range(51)]]

    df = prediction(df, args.output_file)
    plot_rel_rat(dfa, dfb, df, args.output_file)
