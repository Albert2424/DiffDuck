import numpy as np
import pandas as pd
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--results_file', type=str, help='Path to the file containing the DD results. Must be csv', required=True)
parser.add_argument('--output_file', type=str, help='File that will be saved in output directory and will contain the predictions of the DD run.', required=True)
args = parser.parse_args()

def prediction(df,out_file):
    """
    This function takes a pandas dataframe as input and predicts the chain in which the ligand has attached.
    The prediction is based on the confidence scores of each prediction chain.

    Args:
        df (pd.DataFrame): A pandas dataframe containing the prediction results and customer information.
                            The columns of the dataframe must include:
                                'ID', 'Confidence', 'Chain(A=0)(B=1)'

    Returns:
        pd.DataFrame: A pandas dataframe containing the chain ID, prediction, and liability of the prediction.
    """

    pred_df = pd.DataFrame(columns=['ID','Prediction', 'Liability'])
    
    for i in np.unique(df['ID']):

        slice = df[df['ID'] == i].copy()
        slice = slice.reset_index(drop=True)
        # print(slice)
        ci = slice['Confidence']
        c_max = ci[0]
        c_min = ci[len(slice)-1]

        numerator = ci-c_min
        denominator = c_max-c_min

        norm_conf = np.array(numerator)/denominator

        slice['Normal conf'] = norm_conf

        pred_1 = np.sum(np.array(slice['Chain(A=0)(B=1)'])*np.array(slice['Normal conf']))
        pred_1 = pred_1/len(slice) 

        if pred_1 > 0.5:
            pred_df = pd.concat([pred_df if not pred_df.empty else None, pd.DataFrame([{'ID':i,'Prediction':1,'Liability':round(pred_1*100,2)}])], ignore_index=True)
            # print(f'The prediction is chain 1 with certainty of {pred_1*100:.2f}%')
        else:
            pred_df = pd.concat([pred_df if not pred_df.empty else None, pd.DataFrame([{'ID':i,'Prediction':0,'Liability':round((1-pred_1)*100,2)}])], ignore_index=True)
            # print(f'The prediction is chain 0 with certainty of {(1-pred_1)*100:.2f}%')
        # print(pred_1)

        pred_df.to_csv('output/'+out_file)
    return pred_df


def plot_lia_rat(dfa,dfb,pred_df,out):
    import matplotlib.pyplot as plt

    #if DD has failed there will be different values at the dfa and dfb so we must avoid them
    for i in range(len(dfa['Prot ID'])):
        id = f"{dfa['Prot ID'].loc[i]}_{dfb['Prot ID'].loc[i]}_{int(dfa['PubChem CID'].loc[i])}"

        if id in pred_df['ID'].values:
            pass
        else:
            dfa = dfa.drop(i)
            dfb = dfb.drop(i)

    # pred_df = pred_df[pred_df['Prediction'] == 1]
    x = np.array(pred_df['Liability'])
    y = np.log10(np.array(dfa['Ki (nM)'])/np.array(dfb['Ki (nM)']))
    
    plt.figure()
    plt.scatter(x,y)
    plt.xlabel('Liability')
    plt.ylabel(fr'$log(K_i^A/K_i^B)$')
    plt.savefig(f"output/figures/{out.split('.')[0]}.pdf")
    # plt.show()




if __name__ == '__main__':
    print(f'Ranking {args.results_file} values...')
    df = pd.read_csv(args.results_file)
    dfa = pd.read_csv('output/data_A.csv')
    dfb = pd.read_csv('output/data_B.csv')

    df = prediction(df,args.output_file)
    plot_lia_rat(dfa,dfb,df,args.output_file)