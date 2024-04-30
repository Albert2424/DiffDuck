import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from argparse import ArgumentParser
from sklearn import metrics
from sklearn.preprocessing import StandardScaler,PolynomialFeatures,MinMaxScaler
from sklearn.model_selection import cross_val_predict,GridSearchCV
import model_plots as mp
from sklearn import svm
from joblib import dump, load
import datetime
import time
from sklearn.model_selection import ShuffleSplit
import matplotlib.pyplot as plt
import warnings
import pickle
import os

warnings.filterwarnings("ignore")


# Load training data
def load_features(data):
    # Load training data

    train = pd.read_csv(f"{data}/training_data.csv")
    test = pd.read_csv(f"{data}/test_data.csv")

    X_train = train.iloc[:, :-1]
    y_train = train["y_train"]

    X_test = test.iloc[:, :-1]
    y_test = test["y_test"]

    return X_train, y_train, X_test, y_test

def scale(X_train, X_test):
    # Scale the data
    feature_scaler = StandardScaler()
    X_train = feature_scaler.fit_transform(X_train)
    X_test = feature_scaler.transform(X_test)
    return X_train, X_test

def polynomial(X_train, X_test):
    # Make data separable by using polynomial behaviour
    poly = PolynomialFeatures(3)
    X_train = poly.fit_transform(X_train)
    X_test = poly.transform(X_test)
    return X_train, X_test

def single_SVC_ABC(data):

    red = "\033[31m"
    green = "\033[92m"
    end = "\033[0m"


    X_train, y_train, X_test, y_test = load_features(data)
    X_train, X_test = scale(X_train, X_test)

    clf = svm.SVC(C=5, kernel='poly', degree=3,coef0=1).fit(X_train, y_train)

    y_pred = clf.predict(X_test)

    precision_ini = metrics.precision_score(y_test, y_pred, average="macro",labels=np.unique(y_pred))
    recall_ini = metrics.recall_score(y_test, y_pred, average="macro",labels=np.unique(y_pred))
    score = clf.score(X_test,y_test)
    
    print(f"~~~~~~~ {green}SVC metrics{end} ~~~~~~~")
    print(f' Score: {score:.3f}')
    print(f' Precision: {precision_ini:.3f}')
    print(f' Recall: {recall_ini:.3f}')
    print(f" F1: {2/(1/precision_ini+1/recall_ini):.3f}")
    return clf

def grid_SVC_ABC(data,grid_param,model_name,kernel='poly'):

    red = "\033[31m"
    green = "\033[42m"
    end = "\033[0m"

    X_train, y_train, X_test, y_test = load_features(data)
    X_train, X_test = scale(X_train, X_test)

    est = False
    try:
        with open(".grid_params.pkl", "rb") as f:
            loaded_dict = pickle.load(f)

        if loaded_dict == grid_param:
            with open(".et.txt", "r") as file:
                t = file.readlines()
            est = True
        else:
            with open(".grid_params.pkl", "wb") as f:
                pickle.dump(grid_param, f)

    except:
        with open(".grid_params.pkl", "wb") as f:
            pickle.dump(grid_param, f)

    t1 = time.time()
    if est:
        print(f"----> Estimated time: {t} (Current time:{datetime.datetime.now()})")
    else:
        print(f"----> Current time:{datetime.datetime.now()}")

    clf = svm.SVC(kernel=kernel)
    gd_sr = GridSearchCV(estimator=clf,
                        param_grid=grid_param,
                        scoring='balanced_accuracy',
                        cv=5,
                        verbose=1).fit(X_train, y_train)
    t2 = time.time()

    print(f"----> Elapsed time: {str(datetime.timedelta(seconds=round(t2-t1)))}")

    # Collect the best model acording to gd_sr
    best_parameters = gd_sr.best_params_
    best_estimator = gd_sr.best_estimator_

    # Warn the user if the best parameters are the last ones tried (suggest to extend the grid search)
    for i in best_parameters:
        if best_parameters[i] == grid_param[i][-1]:
            try:
                if best_parameters[i] not in [True, False]:
                    float(best_parameters[i])
                    print(
                        f"{red}WARNING{end}: Best {i} is {best_parameters[i]} which is the {red}LAST{end} parameter tried!!!"
                    )
            except ValueError:
                pass

        if best_parameters[i] == grid_param[i][0]:
            try:
                if best_parameters[i] not in [True, False]:
                    float(best_parameters[i])
                    print(
                        f"{red}WARNING{end}: Best {i} is {best_parameters[i]} which is the {red}FIRST{end} parameter tried!!!"
                    )
            except ValueError:
                pass

    best_estimator.fit(X_train, y_train)
    ac = best_estimator.score(X_test, y_test)
    print(f"~~~~~~~ {green}SVC metrics{end} ~~~~~~~")
    print(f' Kernel: {kernel}')
    print(f' Best parameters: {best_parameters}')
    print(f' Accuracy: {ac:.3f}')
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    
    # Save pickle of the model
    dump(
        best_estimator, f"model/{kernel}_{model_name}"
    )  # To load the model: clf = load('model/model.joblib')
    return best_estimator

def single_SVR_logk(data):

    red = "\033[31m"
    green = "\033[92m"
    end = "\033[0m"


    X_train, y_train, X_test, y_test = load_features(data)
    X_train, X_test = scale(X_train, X_test)

    clf = svm.SVR(C=1, kernel='poly', degree=3,coef0=1).fit(X_train, y_train)

    y_pred = clf.predict(X_test)

    # precision_ini = metrics.precision_score(y_test, y_pred, average="macro",labels=np.unique(y_pred))
    # recall_ini = metrics.recall_score(y_test, y_pred, average="macro",labels=np.unique(y_pred))
    score = clf.score(X_test,y_test)
    
    print(f"~~~~~~~ {green}SVR metrics{end} ~~~~~~~")
    print(f' Max Error: {score:.3f}')
    # print(f' Precision: {precision_ini:.3f}')
    # print(f' Recall: {recall_ini:.3f}')
    # print(f" F1: {2/(1/precision_ini+1/recall_ini):.3f}")
    return clf


def grid_SVR_logk(data,grid_param,model_name,kernel='poly'):

    red = "\033[31m"
    green = "\033[42m"
    end = "\033[0m"

    X_train, y_train, X_test, y_test = load_features(data)
    X_train, X_test = scale(X_train, X_test)

    est = False
    try:
        with open(".grid_params.pkl", "rb") as f:
            loaded_dict = pickle.load(f)

        if loaded_dict == grid_param:
            try:
                with open(".et.txt", "r") as file:
                    t = file.readlines()
                est = True
            except:
                pass
        else:
            with open(".grid_params.pkl", "wb") as f:
                pickle.dump(grid_param, f)
            os.remove(".et.txt")

    except:
        with open(".grid_params.pkl", "wb") as f:
            pickle.dump(grid_param, f)

    t1 = time.time()
    if est:
        print(f"----> Estimated time: {t} (Current time:{datetime.datetime.now()})")
    else:
        print(f"----> Current time:{datetime.datetime.now()}")

    clf = svm.SVR(kernel=kernel)
    gd_sr = GridSearchCV(estimator=clf,
                        param_grid=grid_param,
                        scoring='max_error',
                        cv=5,
                        verbose=1).fit(X_train, y_train)
    t2 = time.time()

    print(f"----> Elapsed time: {str(datetime.timedelta(seconds=round(t2-t1)))}")

    # Collect the best model acording to gd_sr
    best_parameters = gd_sr.best_params_
    best_estimator = gd_sr.best_estimator_

    # Warn the user if the best parameters are the last ones tried (suggest to extend the grid search)
    for i in best_parameters:
        if best_parameters[i] == grid_param[i][-1]:
            try:
                if best_parameters[i] not in [True, False]:
                    float(best_parameters[i])
                    print(
                        f"{red}WARNING{end}: Best {i} is {best_parameters[i]} which is the {red}LAST{end} parameter tried!!!"
                    )
            except ValueError:
                pass

        if best_parameters[i] == grid_param[i][0]:
            try:
                if best_parameters[i] not in [True, False]:
                    float(best_parameters[i])
                    print(
                        f"{red}WARNING{end}: Best {i} is {best_parameters[i]} which is the {red}FIRST{end} parameter tried!!!"
                    )
            except ValueError:
                pass

    best_estimator.fit(X_train, y_train)
    ac = best_estimator.score(X_test, y_test)
    print(f"~~~~~~~ {green}SVR metrics{end} ~~~~~~~")
    print(f' Kernel: {kernel}')
    print(f' Best parameters: {best_parameters}')
    print(f' Accuracy: {ac:.3f}')
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    
    # Save pickle of the model
    dump(
        best_estimator, f"model/{kernel}_{model_name}"
    )  # To load the model: clf = load('model/model.joblib')
    return best_estimator


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument(
        "--data",
        type=str,
        help="Path to the directory containing the test and training data.",
        required=True,
    )
    parser.add_argument(
        "--n_samples",
        type=int,
        help="The number of different samples from which the features will be extracted.",
        required=True,
    )
    parser.add_argument(
        "--model_name",
        type=str,
        help="Name of the model that will be saved in a .joblib file.",
        required=True,
    )
    parser.add_argument(
        "--labels",
        type=str,
        help="Kind of labels: options are: 'ABC' or 'logk'. Default is 'logk'",
        default="logk",
        required=False,
    )
    args = parser.parse_args()

    if args.labels == 'ABC':
        clf = single_SVC_ABC(args.data)

        grid_param = {
            'C': [1,2,3,4,5, 6, 7, 8, 9, 10, 11, 12],
            'degree': [2,3,4,5],
            'coef0': [1,2,3,4,5,6,7],
            'max_iter': [500,1000,1500,2000,2500]
        }
        clf = grid_SVC_ABC(args.data,grid_param,args.model_name,kernel='poly')

        grid_param = {
        'C': [0.001,0.01,0.1,0.5,1,2,3,4,5],
        'gamma': [0.001,0.01,0.1,0.5,1,2,3,4,5,6,7,8],
        }

        clf = grid_SVC_ABC(args.data,grid_param,args.model_name,kernel='rbf')
    elif args.labels == 'logk':
        clf = single_SVR_logk(args.data)

        # Polynomial
        print('')
        grid_param = {
            'C': [1,2,3,4,5, 6, 7, 8, 9, 10, 11, 12],
            'degree': [2,3,4,5],
            'coef0': [1,2,3,4,5],
           'max_iter': [1000,1500,2000,2500]
        }
        clf = grid_SVR_logk(args.data,grid_param,args.model_name,kernel='poly')
        print('--> Polynomial grid complete.')

        # Radial Basis Function (RBF)
        print('')
        grid_param = {
        'C': [0.001,0.01,0.1,0.5,1,2,3,4,5],
        'gamma': [0.001,0.01,0.1,0.5,1,2,3,4,5,6,7,8],
        }
        clf = grid_SVR_logk(args.data,grid_param,args.model_name,kernel='rbf')
        print('--> Radial Basis Function (RBF) grid complete.')


    else:
        print(f"Error: {args.labels} is not a valid label. Use 'logk' or 'ABC'")
