import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.cm as cm
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
    "--threshold",
    type=float,
    help="Threshold from which the data has been studied. Default is 0.1",
    default=0.1,
    required=False,
)
parser.add_argument(
    "--predictions",
    type=str,
    help="Path to file containing the DD predictions.",
    required=True,
)
parser.add_argument(
    "--counts",
    type=bool,
    help="Whether to plot or not the counts graph. Default is False",
    required=False,
    default=False,
)
parser.add_argument(
    "--failed_file",
    type=str,
    help="File containing if the DiffDock run seems to fave failed when only running one protein with one ligand. Default is None",
    required=False,
    default=None,
)
args = parser.parse_args()


def plot_ab_count(dfa, dfb, threshold=0.1):
    """
    This function plots the dataframes dfa and dfb, where dfa contains the affinities of proteins A with the same ligand, and dfb contains the affinities of proteins
    B with the same ligand.

    Args:
        dfa (pd.DataFrame): The dataframe containing the affinities of proteins A with the same ligand.
        dfb (pd.DataFrame): The dataframe containing the affinities of proteins B with the same ligand.
        threshold (float, optional): The threshold used for sorting. Defaults to 0.1.
        counts_graph (bool, optional): If True, a histogram of the ligand counts is plotted alongside the error bars. Defaults to False.
        filename (str, optional): The name of the file containing the DD guess. Defaults to None.

    Returns:
        None

    """

    un_lig, counts = np.unique(dfa["SMILES"], return_counts=True)
    col = []

    for lig in dfa["SMILES"]:
        col.append(counts[list(un_lig).index(lig)])

    c_list = cm.viridis(np.linspace(0, 1, 4))

    x = np.array(dfa["Kd (nM)"])
    y = np.array(dfb["Kd (nM)"])

    err_x = np.array(dfa["kd SEM"]) / x
    err_y = np.array(dfb["kd SEM"]) / y

    mask_x = err_x != 0
    mask_y = err_y != 0

    mask = mask_x | mask_y
    err_x = err_x[mask]
    err_y = err_y[mask]

    x_err = x[mask]
    y_err = y[mask]

    x = np.log10(x)
    y = np.log10(y)

    x_err = np.log10(x_err)
    y_err = np.log10(y_err)

    x_values = np.linspace(min(min(x), min(y)) - 5, max(max(x), max(y) + 5), 100)

    plt.figure(figsize=(20, 15))

    plt.errorbar(
        x_err,
        y_err,
        yerr=err_y,
        xerr=err_x,
        alpha=0.3,
        ecolor="r",
        capsize=5,
        fmt="none",
        zorder=0,
    )
    plt.scatter(x, y, c=col, cmap="viridis", s=35)

    plt.plot(x_values, x_values, label="$log(K_i^A) = log(K_i^B)$", color=c_list[0])
    plt.plot(
        x_values,
        x_values + np.log10(threshold),
        color=c_list[1],
        label=f"$log(K_i^B/K_i^A)$ = {threshold}",
    )

    plt.xlabel(rf"$log(K_i^A)$", fontsize=35)
    plt.ylabel(rf"$log(K_i^B)$", fontsize=35)
    plt.xticks(fontsize=35)
    plt.yticks(fontsize=35)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=35)
    cbar.set_label("Ligand counts", fontsize=35)
    plt.legend(fontsize=35)

    plt.xlim([min(min(x), min(x_err)) - 0.2, max(max(x), max(x_err)) + 0.2])
    plt.ylim([min(min(y), min(y_err)) - 0.2, max(max(y), max(y_err)) + 0.2])
    plt.tight_layout()

    plt.savefig("output/figures/ka_kb_count.pdf")
    # plt.show()


def plot_bool_res(dfa, dfb, threshold=0.1, filename=None):
    """
    This function plots the dataframes dfa and dfb, where dfa contains the affinities of proteins A with the same ligand, and dfb contains the affinities of proteins
    B with the same ligand.

    Args:
        dfa (pd.DataFrame): The dataframe containing the affinities of proteins A with the same ligand.
        dfb (pd.DataFrame): The dataframe containing the affinities of proteins B with the same ligand.
        threshold (float, optional): The threshold used for sorting. Defaults to 0.1.
        counts_graph (bool, optional): If True, a histogram of the ligand counts is plotted alongside the error bars. Defaults to False.
        filename (str, optional): The name of the file containing the DD guess. Defaults to None.

    Returns:
        None

    """

    col = pd.read_csv(filename)
    c_list = cm.PiYG_r(np.linspace(0, 1, 2))

    # if DD has failed there will be different values at the dfa and dfb so we must avoid them
    for i in range(len(dfa["Prot ID"])):
        id = f"{dfa['Prot ID'].loc[i]}_{dfb['Prot ID'].loc[i]}_{int(dfa['PubChem CID'].loc[i])}"

        if id in col["ID"].values:
            pass
        else:
            dfa = dfa.drop(i)
            dfb = dfb.drop(i)

    col = np.array(col["Prediction"], dtype=float)
    col[col == 0] = 0.2
    correct = len(col[col == 0.2])
    failed = len(col) - correct

    x = np.array(dfa["Kd (nM)"])
    y = np.array(dfb["Kd (nM)"])

    if len(x) != len(col):
        print(
            f'{filename.split("/")[-1].split(".")[0]}_kakb_pred.pdf will not be generated since the prediction file and the data_A and data_B are not corresponding:{dfa.shape}, {len(col)}'
        )
        return
    err_x = np.array(dfa["kd SEM"]) / x
    err_y = np.array(dfb["kd SEM"]) / y

    mask_x = err_x != 0
    mask_y = err_y != 0

    mask = mask_x | mask_y
    err_x = err_x[mask]
    err_y = err_y[mask]

    x_err = x[mask]
    y_err = y[mask]

    x = np.log10(x)
    y = np.log10(y)

    x_err = np.log10(x_err)
    y_err = np.log10(y_err)

    x_values = np.linspace(min(min(x), min(y)) - 5, max(max(x), max(y) + 5), 100)

    plt.figure(figsize=(20, 15))

    plt.errorbar(
        x_err, y_err, yerr=err_y, xerr=err_x, alpha=0.3, capsize=5, fmt="none", zorder=0
    )
    a = plt.scatter(x, y, c=col, cmap="PiYG_r", s=35, clim=(0, 1))  # 'RdYlGn'

    plt.plot(x_values, x_values, label="$log(K_i^A) = log(K_i^B)$", color=c_list[0])
    plt.plot(
        x_values,
        x_values + np.log10(threshold),
        color=c_list[1],
        label=f"$log(K_i^B/K_i^A)$ = {threshold}",
    )

    plt.xlabel(rf"$log(K_i^A)$", fontsize=35)
    plt.ylabel(rf"$log(K_i^B)$", fontsize=35)
    plt.xticks(fontsize=35)
    plt.yticks(fontsize=35)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=35)

    plt.legend(
        fontsize=35,
        markerscale=4.0,
        handles=a.legend_elements()[0],
        labels=[f"Guessed: {correct}", f"Failed: {failed}"],
    )
    cbar.set_label("DiffDock guess", fontsize=35)
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(["A", "B"])

    plt.xlim([min(min(x), min(x_err)) - 0.2, max(max(x), max(x_err)) + 0.2])
    plt.ylim([min(min(y), min(y_err)) - 0.2, max(max(y), max(y_err)) + 0.2])
    plt.tight_layout()

    plt.savefig(f'output/figures/{filename.split("/")[-1].split(".")[0]}_kakb_pred.pdf')
    # plt.show()


def plot_rel_res(dfa, dfb, threshold=0.1, filename=None):
    """
    This function plots the dataframes dfa and dfb, where dfa contains the affinities of proteins A with the same ligand, and dfb contains the affinities of proteins
    B with the same ligand.

    Args:
        dfa (pd.DataFrame): The dataframe containing the affinities of proteins A with the same ligand.
        dfb (pd.DataFrame): The dataframe containing the affinities of proteins B with the same ligand.
        threshold (float, optional): The threshold used for sorting. Defaults to 0.1.
        counts_graph (bool, optional): If True, a histogram of the ligand counts is plotted alongside the error bars. Defaults to False.
        filename (str, optional): The name of the file containing the DD guess. Defaults to None.

    Returns:
        None

    """
    col = pd.read_csv(filename)
    c_list = cm.cividis_r(np.linspace(0, 1, 2))

    # if DD has failed there will be different values at the dfa and dfb so we must avoid them
    for i in range(len(dfa["Prot ID"])):
        id = f"{dfa['Prot ID'].loc[i]}_{dfb['Prot ID'].loc[i]}_{int(dfa['PubChem CID'].loc[i])}"

        if id in col["ID"].values:
            pass
        else:
            dfa = dfa.drop(i)
            dfb = dfb.drop(i)
    col = np.array(col["Reliability"])
    correct = len(col[col >= 25])
    failed = len(col) - correct

    x = np.array(dfa["Kd (nM)"])
    y = np.array(dfb["Kd (nM)"])

    if len(x) != len(col):
        print(
            f'{filename.split("/")[-1].split(".")[0]}_kakb_rel.pdf will not be generated since the prediction file and the data_A and data_B are not corresponding: {dfa.shape}, {len(col)}'
        )
        return

    err_x = np.array(dfa["kd SEM"]) / x
    err_y = np.array(dfb["kd SEM"]) / y

    mask_x = err_x != 0
    mask_y = err_y != 0

    mask = mask_x | mask_y
    err_x = err_x[mask]
    err_y = err_y[mask]

    x_err = x[mask]
    y_err = y[mask]

    x = np.log10(x)
    y = np.log10(y)

    x_err = np.log10(x_err)
    y_err = np.log10(y_err)

    x_values = np.linspace(min(min(x), min(y)) - 5, max(max(x), max(y) + 5), 100)

    plt.figure(figsize=(20, 15))

    plt.errorbar(
        x_err, y_err, yerr=err_y, xerr=err_x, alpha=0.3, capsize=5, fmt="none", zorder=0
    )
    a = plt.scatter(x, y, c=col, cmap="cividis_r", s=35, clim=(0, 100))  # 'RdYlGn'

    plt.plot(x_values, x_values, label="$log(K_i^A) = log(K_i^B)$", color=c_list[0])
    plt.plot(
        x_values,
        x_values + np.log10(threshold),
        color=c_list[1],
        label=f"$log(K_i^B/K_i^A)$ = {threshold}",
    )

    plt.xlabel(rf"$log(K_i^A)$", fontsize=35)
    plt.ylabel(rf"$log(K_i^B)$", fontsize=35)
    plt.xticks(fontsize=35)
    plt.yticks(fontsize=35)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=35)

    plt.legend(
        fontsize=35,
        markerscale=4.0,
        handles=[a.legend_elements()[0][-1], a.legend_elements()[0][0]],
        labels=[f"Reliability > 25%: {correct}", f"Reliability < 25%: {failed}"],
    )
    cbar.set_label("DiffDock guess reliability", fontsize=35)
    cbar.set_ticks([50, 75, 100])
    cbar.set_ticklabels(["50%", "75%", "100%"])

    plt.xlim([min(min(x), min(x_err)) - 0.2, max(max(x), max(x_err)) + 0.2])
    plt.ylim([min(min(y), min(y_err)) - 0.2, max(max(y), max(y_err)) + 0.2])
    plt.tight_layout()

    plt.savefig(f'output/figures/{filename.split("/")[-1].split(".")[0]}_kakb_rel.pdf')
    # plt.show()


def plot_failed(dfa, dfb, failed, threshold=0.1, counts=False, predictions=None):

    try:
        fail_df = pd.read_csv(failed)
    except:
        print(f"Wrong file was inputed: {failed}")
        return

    pred = pd.read_csv(predictions)
    pred = pred[["ID", "Prediction", "Reliability"]]

    # if DD has failed there will be different values at the dfa and dfb so we must avoid them
    order = []
    for i in range(len(dfa["Prot ID"])):
        id = f"{dfa['Prot ID'].loc[i]}_{dfb['Prot ID'].loc[i]}_{int(dfa['PubChem CID'].loc[i])}"
        order.append(id)
        if id in pred["ID"].values:
            pass
        else:
            dfa = dfa.drop(i)
            dfb = dfb.drop(i)

    pred = pred.sort_values(by=["ID"])

    # Sort the predictions to match the dfa and dfb
    pred["ID"] = pd.Categorical(pred["ID"], categories=order, ordered=True)
    pred = pred.sort_values(by="ID")

    pred = pred.reset_index()

    for i in range(len(fail_df["Fail"])):
        f = fail_df.loc[i, "Fail"]
        prot = fail_df.loc[i, "Protein"]
        cid = fail_df.loc[i, "PubChem CID"]

        if f == 1:
            mask_a = dfa["Prot ID"].eq(prot) & dfa["PubChem CID"].eq(cid)
            idx_to_drop = dfa.index[mask_a]
            dfa = dfa.drop(idx_to_drop)
            dfb = dfb.drop(idx_to_drop)
            pred = pred.drop(idx_to_drop)

            mask_b = dfb["Prot ID"].eq(prot) & dfb["PubChem CID"].eq(cid)
            idx_to_drop = dfa.index[mask_b]
            dfa = dfa.drop(idx_to_drop)
            dfb = dfb.drop(idx_to_drop)
            pred = pred.drop(idx_to_drop)

            dfa = dfa.reset_index(drop=True)
            dfb = dfb.reset_index(drop=True)
            pred = pred.reset_index(drop=True)

    pred.to_csv(predictions.split(".")[0] + "_no_fails.csv")

    predictions = predictions.split(".")[0] + "_no_fails.csv"
    if counts:
        plot_ab_count(dfa, dfb, threshold=threshold)
    plot_bool_res(dfa, dfb, threshold=threshold, filename=predictions)
    plot_rel_res(dfa, dfb, threshold=threshold, filename=predictions)


if __name__ == "__main__":
    print(f"Plotting {args.predictions}")
    dfa = pd.read_csv(args.data_a)
    dfb = pd.read_csv(args.data_b)

    if args.failed_file != "None":
        print(
            f"Using {args.failed_file} to remove pairs that have individually failed with their ligand..."
        )
        print("")
        plot_failed(
            dfa,
            dfb,
            args.failed_file,
            threshold=args.threshold,
            counts=args.counts,
            predictions=args.predictions,
        )

    else:
        if args.counts:
            plot_ab_count(dfa, dfb, threshold=args.threshold)
        plot_bool_res(dfa, dfb, threshold=args.threshold, filename=args.predictions)
        plot_rel_res(dfa, dfb, threshold=args.threshold, filename=args.predictions)
