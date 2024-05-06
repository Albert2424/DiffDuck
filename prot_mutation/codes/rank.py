import numpy as np
import pandas as pd
from argparse import ArgumentParser
import os

def args_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "--run_name",
        type=str,
        help="Name of the job.",
        required=True,
    )
    args = parser.parse_args()
    return args

def rank_docking_all(run_name):
    results = [f'result_{run_name}_AF.csv',f'result_{run_name}_DF.csv',f'result_{run_name}_OF.csv']
    rank_total = np.array([])
    fold = []
    for r in results:
        rank = []
        res = pd.read_csv(f'results/{r}')

        # samples = len(res[res['ID'] == res['ID'].iloc[0]])
        
        ids = np.unique(res['ID'].values)
        if len(rank_total) == 0:
            rank_total = np.empty(len(ids))

        
        for id in ids:
            rank.append(np.sum(res[res['ID'] == id]['Chain(A=0)(B=1)'].values[:5])/5)
            fold.append(r.split('.')[0].split('_')[-1])
        # for i in range(len(rank)):
        #     print(ids[i],rank[i])

        rank_total += np.array(rank) 

    ranked = pd.DataFrame()


    ranked['ID'] = ids
    ranked['Rank'] = np.round(rank_total/3,4)
    ranked['Folding'] = ['all' for i in range(len(ids))]

    ranked.to_csv(f'results/ranked_{run_name}_all.csv', index=False)




def rank_one(run_name, fold='AF'):
    results = f'results/result_{run_name}_{fold}.csv'

    res = pd.read_csv(results)
    # samples = len(res[res['ID'] == res['ID'].iloc[0]])
    
    ids = np.unique(res['ID'].values)

    rank = []
    for id in ids:
        rank.append(np.sum(res[res['ID'] == id]['Chain(A=0)(B=1)'].values[:5])/5)

    # for i in range(len(rank)):
    #     print(ids[i],rank[i])

    ranked = pd.DataFrame()

    ranked['ID'] = ids
    ranked['Rank'] = rank
    ranked['Folding'] = [fold for i in range(len(ids))]

    ranked.to_csv(f'results/ranked_{run_name}_{fold}.csv', index=False)


if __name__ == '__main__':
    args = args_parser()
    print(f" Ranking docking results from {args.run_name}...")
    rank_one(args.run_name)
    rank_one(args.run_name, fold='DF')
    rank_one(args.run_name, fold='OF')
    rank_docking_all(args.run_name)

