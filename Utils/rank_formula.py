import numpy as np
import pandas as pd

def prediction(df):
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
            pred_df = pd.concat([pred_df, pd.DataFrame([{'ID':i,'Prediction':1,'Liability':round(pred_1*100,2)}])], ignore_index=True)
            print(f'The prediction is chain 1 with certainty of {pred_1*100:.2f}%')
        else:
            pred_df = pd.concat([pred_df, pd.DataFrame([{'ID':i,'Prediction':0,'Liability':round((1-pred_1)*100,2)}])], ignore_index=True)
            print(f'The prediction is chain 0 with certainty of {(1-pred_1)*100:.2f}%')
        # print(pred_1)

    return pred_df


def plot_lia_rat(dfa,dfb,pred_df):
    import matplotlib.pyplot as plt

    x = np.array(pred_df['Liability'])
    y = np.log(np.array(dfa['Ki (nM)'])/np.array(dfb['Ki (nM)']))
    
    plt.figure()
    plt.scatter(x,y)
    plt.xlabel('Liability')
    plt.ylabel(fr'$log(K_i^A/K_i^B)$')
    plt.show()



if __name__ == '__main__':
    df = pd.read_csv('results.csv')
    dfa = pd.read_csv('data_A.csv')
    dfb = pd.read_csv('data_B.csv')

    df = df[df['Confidence'] != 0]

    df = df.reset_index(drop=True)

    df = prediction(df)
    # plot_lia_rat(dfa,dfb,df)