import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.cm as cm

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
        c_list = cm.PiYG_r(np.linspace(0, 1, 2))

         #if DD has failed there will be different values at the dfa and dfb so we must avoid them
        for i in range(len(dfa['Prot ID'])):
            id = f"{dfa['Prot ID'].loc[i]}_{dfb['Prot ID'].loc[i]}_{int(dfa['PubChem CID'].loc[i])}"

            if id in col['ID'].values:
                pass
            else:
                dfa = dfa.drop(i)
                dfb = dfb.drop(i)

        col = list(col['Prediction'])
        
        x = np.array(dfa['Ki (nM)'])
        y = np.array(dfb['Ki (nM)'])

    err_x = np.array(dfa['ki SEM'])/x
    err_y = np.array(dfb['ki SEM'])/y

    mask_x = err_x != 0
    mask_y = err_y!= 0

    mask = mask_x | mask_y
    err_x = err_x[mask]
    err_y = err_y[mask]

    x_err = x[mask]
    y_err = y[mask]

    x = np.log(x)
    y = np.log(y)

    x_err = np.log(x_err)
    y_err = np.log(y_err)

    x_values = np.linspace(min(min(x),min(y)), max(max(x),max(y)), 100)

    plt.figure(figsize=(20,15))
    if counts_graph:
        plt.errorbar(x_err, y_err, yerr=err_y, xerr=err_x, ecolor='r', capsize=5,fmt='none', zorder=0)
        plt.scatter(x,y, c=col,cmap='viridis')
    else:
        plt.errorbar(x_err, y_err, yerr=err_y, xerr=err_x, ecolor='r', capsize=5,fmt='none', zorder=0)
        plt.scatter(x,y, c=col,cmap='PiYG_r') # 'RdYlGn'

    plt.plot(x_values, x_values, label='$log(K_i^A) = log(K_i^B)$', color=c_list[0])
    plt.plot(x_values, x_values+np.log(threshold), color=c_list[1] , label=f'$log(K_i^B/K_i^A)$ = {threshold}')
    
    plt.xlabel(fr'$log(K_i^A)$', fontsize=35)
    plt.ylabel(fr'$log(K_i^B)$', fontsize=35)
    plt.legend(fontsize=35)
    plt.xticks(fontsize=35)
    plt.yticks(fontsize=35)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=35)
    if counts_graph:
        cbar.set_label('Ligand counts', fontsize=35)
    else:
        cbar.set_label('DiffDock guess', fontsize=35)
        cbar.set_ticks([0,1])
        cbar.set_ticklabels(['A','B'])
    plt.tight_layout()
    if counts_graph:
        plt.savefig('output/figures/ka_kb_count.pdf')
    else:
        plt.savefig('output/figures/ka_kb_pred.pdf')
    plt.show()


dfa = pd.read_csv('output/data_A.csv')
dfb = pd.read_csv('output/data_B.csv')
thresh = 0.1

plot_ab(dfa,dfb, counts_graph=True,threshold=thresh)
plot_ab(dfa,dfb, counts_graph=False,threshold=thresh,filename='output/pred_results_DD_DF.csv')
