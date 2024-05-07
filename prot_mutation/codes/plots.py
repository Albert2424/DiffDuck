import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.cm as cm
from argparse import ArgumentParser
import seaborn as sns
import os
from pandas.api.types import CategoricalDtype

def args_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "--run_name",
        type=str,
        help="Name of the job.",
        required=True,
    )
    parser.add_argument(
        "--AF",
        type=int,
        help="whether AF was used or not. Default is 1 (true)",
        required=False,
        default=1
    )
    parser.add_argument(
        "--DF",
        type=int,
        help="whether DF was used or not. Default is 1 (true)",
        required=False,
        default=1
    )
    parser.add_argument(
        "--OF",
        type=int,
        help="whether OF was used or not. Default is 1 (true)",
        required=False,
        default=1
    )
    args = parser.parse_args()
    return args

def plot_best_chain(run_name, fold='all'):

    df = pd.read_csv(f"results/ranked_{run_name}_{fold}.csv")
    ids = df['ID'].values
    rank = df['Rank'].values

    prot1 = []
    prot2 = []

    for i in ids:
        id = i.split('_')
        prot1.append(id[0])
        prot2.append(id[1])

    df['Prot 1'] = prot1
    df['Prot 2'] = prot2

    c_list = cm.cividis(np.linspace(0, 1, 2))

    # cat1 = CategoricalDtype(categories=prot1, ordered=True)
    # df["ID"] = df["ID"].astype(cat1)

    data = df.pivot(index="Prot 2", columns="Prot 1", values="Rank")
    data = data.iloc[data.isnull().sum(axis=1).mul(-1).argsort()] # make sure the data is sorted

    # data = data.reindex(data.sum().sort_values(ascending=False).index, axis=1)
 
    plt.figure()
    hm = sns.heatmap(
        data,
        cmap="cividis",
        cbar_kws={"label": "Rank"},
        linewidth=.5,
        annot=True,
        vmin=0,
        vmax=1,
    )
    bbox = dict(boxstyle="round", ec='black', fc=c_list[0], alpha=0.7)
    plt.setp(hm.get_xticklabels(), bbox=bbox)
    bbox = dict(boxstyle="round", ec='black', fc=c_list[1], alpha=0.7)
    plt.setp(hm.get_yticklabels(), bbox=bbox)

    plt.title(f'{run_name} ({fold})')

    plt.tight_layout()
    plt.savefig(f'results/figures/rank_{run_name}_{fold}.png')
    # plt.show()


if __name__ == "__main__":
    args = args_parser()
    print(f'Plotting ranked predictions from {args.run_name}')

    if args.AF == 1:
        plot_best_chain(args.run_name, 'AF')
    if args.DF == 1:
        plot_best_chain(args.run_name, 'DF')
    if args.OF == 1:
        plot_best_chain(args.run_name, 'OF')

    plot_best_chain(args.run_name, 'all')
    # plot_best_chain('FABP', 'all')
