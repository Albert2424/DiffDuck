import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.cm as cm
from argparse import ArgumentParser



def plot_bool_res(dfa, dfb, threshold=0.1, predictions=None, plot_th=0):
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

    # Load the predictions once the DD results have been ranked.
    pred = pd.read_csv(predictions)

    # color list for ploting lines
    c_list = cm.PiYG_r(np.linspace(0, 1, 2))

    # Select from data_A and data_B only the rows that are in the validation set.
    dfa_aux = pd.DataFrame()
    dfb_aux = pd.DataFrame()

    for i, id in enumerate(pred["Prot A ID"]):

        # Select the row that corresponds to prot A and prot B with the ligand
        sa = dfa[dfa["Prot ID"] == id]
        sb = dfb[dfa["Prot ID"] == id]

        sa = sa[sa["PubChem CID"] == pred["Ligand ID"][i]]
        sb = sb[sb["PubChem CID"] == pred["Ligand ID"][i]]

        sa = sa[sb["Prot ID"] == pred["Prot B ID"][i]]
        sb = sb[sb["Prot ID"] == pred["Prot B ID"][i]]

        # Append the row in the auxiliary dataframe
        dfa_aux = pd.concat(
            [
                dfa_aux if not dfa_aux.empty else None,
                sa,
            ],
            ignore_index=True,
        )
        dfb_aux = pd.concat(
            [
                dfb_aux if not dfb_aux.empty else None,
                sb,
            ],
            ignore_index=True,
        )

    # Now dfa and dfb contain only the selected rows.
    dfa = dfa_aux
    dfb = dfb_aux

    # In case that the threshold is different than 1 (no th)
    if threshold != 1:
        mask = dfa[dfa["Kd (nM)"] / dfb["Kd (nM)"] >= threshold]
        pred_th = pred.drop(mask.index)  # select the points below the th

        truth_th = pred_th["Experimental Pred"]
        pred_th = np.array(pred_th["Prediction"], dtype=float)
        pred_th[pred_th == 0] = 0.2
        mask = np.where(pred_th == truth_th)

        # correct and failed below the th
        correct_th = len(pred_th[mask])
        failed_th = len(pred_th) - correct_th

        # ambiguous values that have been correctly labeled as ambiguous
        amb_th = 0
        for i, val in enumerate(pred_th):
            if val == 0.7:
                if i in mask[0]:
                    amb_th += 1
                else:
                    pred_th[i] = 1

        # only plot the points below the th.
        if plot_th == 1:
            dfa = dfa.drop(mask.index)
            dfb = dfb.drop(mask.index)
            pred = pred.drop(mask.index)

    # obtain the correctly and incorrectly predicted labels
    # 0.2 for correct, 1.0 for failed and 0.7 for ambiguous
    truth = pred["Experimental Pred"]
    pred = np.array(pred["Prediction"], dtype=float)
    pred[pred == 0] = 0.2
    truth[truth == 0] = 0.2
    # correct = len(pred[pred == 0.2])
    truth_mask = np.where(pred == truth)
    correct = len(pred[truth_mask])
    failed = len(pred) - correct

    # print(correct/(correct+failed))

    # select X and Y for plotting
    x = np.array(dfb["Kd (nM)"])
    y = np.array(dfa["Kd (nM)"])

    # split into different sets (A, B, ambiguous)
    x_amb = x[pred == 0.7]
    x_a = x[pred == 0.2]
    x_b = x[pred == 1]

    y_amb = y[pred == 0.7]
    y_a = y[pred == 0.2]
    y_b = y[pred == 1]

    mask_amb = pred == 0.7
    mask_a = pred == 0.2
    mask_b = pred == 1

    # Relabel the correct and incorrect values of every split
    amb_c = 0
    amb_f = 0
    a_c = 0
    a_f = 0
    b_c = 0
    b_f = 0
    for i, val in enumerate(pred):
        if val == 0.7:
            if i in truth_mask[0]:
                amb_c += 1
                pred[i] = 0.2
                # print(dfa["Kd (nM)"].values[i]/dfb["Kd (nM)"].values[i])
            else:
                pred[i] = 1
                amb_f += 1
                # print(dfa["Kd (nM)"].values[i]/dfb["Kd (nM)"].values[i])

        elif val == 0.2:
            if i not in truth_mask[0]:
                pred[i] = 1
                a_f += 1
            else:
                pred[i] = 0.2
                a_c += 1
        else:
            if i not in truth_mask[0]:
                pred[i] = 1
                b_f += 1
            else:
                pred[i] = 0.2
                b_c += 1

    red = '\033[31m'
    green = '\033[92m'
    end = '\033[0m'
    print(f'{green}Correct:{end}\n\tTrue A: {a_c}\n\tTrue B: {b_c}\n\tTrue Ambiguous: {amb_c}\n\n{red}Failed:{end}\n\tFalse A: {a_f}\n\tFalse B: {b_f}\n\tFalse Ambiguous: {amb_f}')
    print('---------------------------------------------------\n')

    # select the colouring of every set
    pred_amb = pred[mask_amb]
    pred_a = pred[mask_a]
    pred_b = pred[mask_b]

    # Make sure that the prediction file and the data file are from the same run
    if len(x) != len(pred):
        print(
            f'{predictions.split("/")[-1].split(".")[0]}_kakb_pred.png will not be generated since the prediction file and the data_A and data_B are not corresponding:{dfa.shape}, {len(pred)}'
        )
        return

    # error for every point
    err_x = np.array(dfb["kd SEM"]) / x
    err_y = np.array(dfa["kd SEM"]) / y

    mask_x = err_x != 0
    mask_y = err_y != 0

    mask = mask_x | mask_y
    err_x = err_x[mask]
    err_y = err_y[mask]

    x_err = x[mask]
    y_err = y[mask]

    x = np.log10(x)
    y = np.log10(y)

    x_amb = np.log10(x_amb)
    x_a = np.log10(x_a)
    x_b = np.log10(x_b)

    y_amb = np.log10(y_amb)
    y_a = np.log10(y_a)
    y_b = np.log10(y_b)

    x_err = np.log10(x_err)
    y_err = np.log10(y_err)

    x_values = np.linspace(min(min(x), min(y)) - 5, max(max(x), max(y) + 5), 100)

    plt.figure(figsize=(20, 15))

    plt.errorbar(
        x_err, y_err, yerr=err_y, xerr=err_x, alpha=0.3, capsize=5, fmt="none", zorder=0
    )
    print(f"Correct: {correct}, failed: {failed}")
    s1 = plt.scatter(
        x_a,
        y_a,
        c=pred_a,
        cmap="PiYG_r",
        s=35,
        clim=(0, 1),
        marker="D",
        label="Predicted as A",
    )
    s2 = plt.scatter(
        x_b,
        y_b,
        c=pred_b,
        cmap="PiYG_r",
        s=35,
        clim=(0, 1),
        marker="X",
        label="Predicted as B",
    )
    s3 = plt.scatter(
        x_amb,
        y_amb,
        c=pred_amb,
        cmap="PiYG_r",
        s=35,
        clim=(0, 1),
        label="Predicted as Ambiguous",
    )

    leg1 = plt.legend(handles=[s1,s2,s3],labels=["Predicted as A","Predicted as B","Predicted as Ambiguous"],fontsize=35, markerscale=4.0)
    # a = plt.scatter(x, y, c=pred, cmap="PiYG_r", s=35, clim=(0, 1))  # 'RdYlGn'
    plt.gca().add_artist(leg1)

    l1, = plt.plot(
        x_values, x_values, label="$log(K_D^A) = log(K_D^B)$", color=c_list[0]
    )
    l2, = plt.plot(
        x_values,
        x_values + np.log10(threshold),
        color=c_list[1],
        label=f"$log(K_D^A/K_D^B) = log({threshold}$)",
    )

    leg2 = plt.legend(handles=[l1,l2],labels=["$log(K_D^A) = log(K_D^B)$",f"$log(K_D^A/K_D^B) = log({threshold}$)"],loc='upper center', bbox_to_anchor=(0.5, -0.1),
          fancybox=True, shadow=True, ncol=2, fontsize=35)
    # plt.gca().add_artist(leg2)
    plt.xlabel(rf"$log(K_D^B)$", fontsize=35)
    plt.ylabel(rf"$log(K_D^A)$", fontsize=35)
    plt.xticks(fontsize=35)
    plt.yticks(fontsize=35)

    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=35)
    cbar.set_label("DiffDock guess", fontsize=35)
    cbar.set_ticks([0.2, 0.95])
    cbar.set_ticklabels(["Correct", "Failed"], rotation=90)

    try:
        min_value = min(min(min(x), min(x_err)), min(min(y), min(y_err))) - 0.2
        max_value = max(max(max(x), max(x_err)), max(max(y), max(y_err))) + 0.2
    except ValueError:
        min_value = min(min(x), min(y)) - 0.2
        max_value = max(max(x), max(y)) + 0.2

    plt.xlim([min_value, max_value])
    plt.ylim([min_value, max_value])

    plt.tight_layout()
    plt.savefig(
        f'output/figures/{predictions.split("/")[-1].split(".")[0]}_kakb_pred.png'
    )
    # plt.show()


def plot_logk_pred(pred):

    pred_df = pd.read_csv(pred)

    y_true = np.array(pred_df["Prediction"], dtype=float)
    y_pred = np.array(pred_df["Experimental Pred"], dtype=float)

    c_list = cm.PiYG_r(np.linspace(0.2, 1, 2))
    threshold = 0.1
    x_values = np.linspace(min(min(y_pred), min(y_true)) - 5, max(max(y_pred), max(y_true) + 5), 100)

    plt.figure(figsize=(20, 15))
    plt.scatter(y_pred,y_true)

    plt.xlabel(r"$log(\frac{K_D^A}{K_D^B})$ Predicted", fontsize=35)
    plt.ylabel(r"$log(\frac{K_D^A}{K_D^B})$ Real", fontsize=35)
    plt.xticks(fontsize=35)
    plt.yticks(fontsize=35)

    l1, = plt.plot(
        x_values, x_values, label="Prediction = True value", color=c_list[0]
    )
    l2, = plt.plot(
        x_values,
        x_values + np.log10(threshold),
        color=c_list[1],
        label=f"$log(K_D^A/K_D^B) = log({threshold}$)",
        linestyle="--", 
    )
    l3, = plt.plot(
        x_values,
        x_values - np.log10(threshold),
        color=c_list[1],
        label=f"$log(K_D^A/K_D^B) = log({threshold}$)",
        linestyle="--",
    )


    leg2 = plt.legend(handles=[l1,l2],labels=["$log(K_D^A) = log(K_D^B)$",f"$log(K_D^A/K_D^B) = log({threshold}$)"],loc='upper center', bbox_to_anchor=(0.5, -0.1),
          fancybox=True, shadow=True, ncol=2, fontsize=35)
    

    min_value = min(min(y_pred), min(y_true)) - 0.2
    max_value = max(max(y_pred), max(y_true)) + 0.2


    plt.fill_between(x_values, x_values + np.log10(threshold), x_values - np.log10(threshold), color=c_list[0], alpha=0.05)

    plt.xlim([min_value, max_value])
    plt.ylim([min_value, max_value])

    plt.tight_layout()
    plt.savefig(
        f'output/figures/{pred.split("/")[-1].split(".")[0]}_logk.png'
    )
    

if __name__ == "__main__":
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
        "--failed_file",
        type=str,
        help="File containing if the DiffDock run seems to fave failed when only running one protein with one ligand. Default is None",
        required=False,
        default=None,
    )
    parser.add_argument(
        "--plot_th",
        type=int,
        help="Whether or not to only plot the points below the threshold. Default is False.",
        required=False,
        default=0,
    )
    parser.add_argument(
        "--labels",
        type=str,
        help="Kind of labels: options are: 'ABC' or 'logk'. Default is 'logk'",
        default="logk",
        required=False,
    )
    args = parser.parse_args()

    print(f"Plotting {args.predictions}")
    dfa = pd.read_csv(args.data_a)
    dfb = pd.read_csv(args.data_b)
    
    if args.labels == 'ABC':
        plot_bool_res(
            dfa,
            dfb,
            threshold=args.threshold,
            predictions=args.predictions,
            plot_th=args.plot_th,
        )
    elif args.labels == 'logk':
        plot_logk_pred(args.predictions)
