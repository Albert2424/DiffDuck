import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.cm as cm
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--data_a', type=str, help='Path to the file containing the proteins and ligands with their affinities of proteins A (more affine). Must be csv', required=True)
parser.add_argument('--data_b', type=str, help='Path to the file containing the proteins and ligands with their affinities of proteins B (less affine). Must be csv', required=True)
parser.add_argument('--threshold', type=float, help='Threshold from which the data has been studied. Default is 0.1', default=0.1, required=False)
parser.add_argument('--predictions', type=str, help='Path to file containing the DD predictions.',required=True)
args = parser.parse_args()

def plot_ab(dfa,dfb, threshold=0.1, counts_graph=False, filename=None):

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

    

    if counts_graph:
        un_lig, counts = np.unique(dfa['SMILES'],return_counts=True)
        col = []

        for lig in dfa['SMILES']:
            col.append(counts[list(un_lig).index(lig)])
        
        c_list = cm.viridis(np.linspace(0, 1, 4))

        x = np.array(dfa['Ki (nM)'])
        y = np.array(dfb['Ki (nM)'])

    else:
        col = pd.read_csv(filename)
        c_list = cm.PiYG_r(np.linspace(0, 1, 4))

         #if DD has failed there will be different values at the dfa and dfb so we must avoid them
        for i in range(len(dfa['Prot ID'])):
            id = f"{dfa['Prot ID'].loc[i]}_{dfb['Prot ID'].loc[i]}_{int(dfa['PubChem CID'].loc[i])}"

            if id in col['ID'].values:
                pass
            else:
                dfa = dfa.drop(i)
                dfb = dfb.drop(i)

        col = np.array(col['Prediction'])

        correct = len(col[col==0])
        failed = len(col)-correct
        
        x = np.array(dfa['Ki (nM)'])
        y = np.array(dfb['Ki (nM)'])

    err_x = np.array(dfa['ki SEM'])/x
    err_y = np.array(dfb['ki SEM'])/y

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

    x_values = np.linspace(min(min(x),min(y))-5, max(max(x),max(y)+5), 100)

    plt.figure(figsize=(20,15))
    if counts_graph:
        plt.errorbar(x_err, y_err, yerr=err_y, xerr=err_x, alpha=0.3, ecolor='r', capsize=5,fmt='none', zorder=0)
        plt.scatter(x,y, c=col,cmap='viridis', s=35)
    else:
        plt.errorbar(x_err, y_err, yerr=err_y, xerr=err_x, alpha=0.3, capsize=5,fmt='none', zorder=0)
        a = plt.scatter(x,y, c=col,cmap='PiYG_r', s=35) # 'RdYlGn'

    plt.plot(x_values, x_values, label='$log(K_i^A) = log(K_i^B)$', color=c_list[0])
    plt.plot(x_values, x_values+np.log10(threshold), color=c_list[1] , label=f'$log(K_i^B/K_i^A)$ = {threshold}')
    
    plt.xlabel(fr'$log(K_i^A)$', fontsize=35)
    plt.ylabel(fr'$log(K_i^B)$', fontsize=35)
    plt.xticks(fontsize=35)
    plt.yticks(fontsize=35)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=35)
    if counts_graph:
        cbar.set_label('Ligand counts', fontsize=35)
        plt.legend(fontsize=35)
    else:
        plt.legend(fontsize=35,markerscale=4.,handles=a.legend_elements()[0], labels=[f'Guessed: {correct}',f'Failed: {failed}'])
        cbar.set_label('DiffDock guess', fontsize=35)
        cbar.set_ticks([0,1])
        cbar.set_ticklabels(['A','B'])

    plt.xlim([min(min(x),min(x_err))-0.2,max(max(x),max(x_err))+0.2])
    plt.ylim([min(min(y),min(y_err))-0.2,max(max(y),max(y_err))+0.2])
    plt.tight_layout()
    if counts_graph:
        plt.savefig('output/figures/ka_kb_count.pdf')
    else:
        plt.savefig(f'output/figures/{filename.split("/")[-1].split(".")[0]}_pred.pdf')
    # plt.show()


if __name__ == '__main__':
    print(f'Plotting {args.predictions}')
    dfa = pd.read_csv(args.data_a)
    dfb = pd.read_csv(args.data_b)
    thresh = .1

    plot_ab(dfa,dfb, counts_graph=True,threshold=thresh)
    plot_ab(dfa,dfb, counts_graph=False,threshold=thresh,filename=args.predictions)
