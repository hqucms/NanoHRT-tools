from __future__ import print_function

import os
import shutil

import xgboost as xgb

import glob
import uproot
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import argparse
parser = argparse.ArgumentParser('xgboost training and application')
parser.add_argument('--train', action='store_true', help='Do training. Otherwise do testing.')
parser.add_argument('--gpu', default='', help='gpu number')
parser.add_argument('--xml-output', default=None, help='Convert xgboost model to TMVA XML.')
parser.add_argument('-i', '--input', help='Input file')
parser.add_argument('-o', '--output', help='Output file')
args = parser.parse_args()

# ---------------------------------------------
train_val_files = glob.glob('/eos/cms/store/user/hqu/resTop/20220711_resTop/pieces/TTToSemiLeptonic_*_tree.root')
train_val_files = train_val_files[:2]  # FIXME: use only 2 files for a quick testing
train_val_split = 0.5

model_file = 'models/xgb_train_resTop.model'

wgtvar = None
label_var = 'isHadTop'

train_vars = (
    # j1
    'j1_pt',
    'j1_mass',
    'j1_deepFlavB',
    'j1_deepFlavCvL',
    'j1_deepFlavCvB',
    'j1_deepFlavQG',

    # j2
    'j2_pt',
    'j2_mass',
    'j2_deepFlavB',
    'j2_deepFlavCvL',
    'j2_deepFlavCvB',
    'j2_deepFlavQG',

    # j3
    'j3_pt',
    'j3_mass',
    'j3_deepFlavB',
    'j3_deepFlavCvL',
    'j3_deepFlavCvB',
    'j3_deepFlavQG',

    # di-jet
    'j12_deltaR',
    'j12_mass',
    'j13_deltaR',
    'j13_mass',
    'j23_deltaR',
    'j23_mass',

    # tri-jet
    'topcand_mass',
)


obs_vars = (
)

all_vars = set(train_vars + obs_vars)

# ---------------------------------------------


def make_dmatrix(filepaths, predict=False, train_val_split=0.5):
    if predict:
        raise NotImplemented("...predict...")
    else:
        dfs = []
        for f in filepaths:
            tree = uproot.open(f)['Friends']
            if uproot.__version__[0] == '3':
                df = tree.pandas.df(branches=list(train_vars) + [label_var])
            else:
                df = tree.arrays(filter_name=list(train_vars) + [label_var], library='pd')
            dfs.append(df)
        a = pd.concat(dfs)
        n_train = int(a.shape[0] * train_val_split)
        X = np.stack([a[var] for var in train_vars], axis=1)
        y = a[label_var]
        assert X.shape[0] == y.shape[0]
        d_train = xgb.DMatrix(X[:n_train], y[:n_train], feature_names=train_vars)
        d_val = xgb.DMatrix(X[n_train:], y[n_train:], feature_names=train_vars)
        return d_train, d_val


def plotROC(y_score, y_true, sample_weight=None):
    from sklearn.metrics import auc, roc_curve
    fpr, tpr, _ = roc_curve(y_true, y_score, sample_weight=sample_weight)
    roc_auc = auc(fpr, tpr)

    plt.figure()
    legend = '%s (area = %0.4f)' % ('res-top', roc_auc)
    print(legend)
    plt.plot(tpr, 1 - fpr, label=legend)
    plt.plot([0, 1], [1, 0], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0, 1])
    plt.xlabel('Signal Efficiency')
    plt.ylabel('Background rejection')
    plt.legend(loc='best')
#     plt.yscale('log')
    plt.grid()


if args.xml_output:
    from xml_helper import toXML
    toXML(train_vars, model_file, args.xml_output)
    import sys
    sys.exit()

if args.train:
    model_dir = os.path.dirname(model_file)
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
    shutil.copy(__file__, model_dir)

    print('Start loading files...')
    d_train, d_test = make_dmatrix(train_val_files, train_val_split=train_val_split)

    print('Using training files with %d events:\n%s\n...' % (d_train.num_row(), train_val_files))
    print('Using validation files with %d events:\n%s\n...' % (d_test.num_row(), train_val_files))

    print('Training vars:\n' + '\n'.join(train_vars))
    print('%d vars in total' % len(train_vars))

    # # setup parameters for xgboost
    param = {}
    param['objective'] = 'binary:logistic'
    param['eval_metric'] = ['error', 'auc', 'logloss']
    param['min_child_weight'] = 1e-4
    # param['gamma']=0.01
    param['eta'] = 0.1
    param['max_depth'] = 6
    # param['colsample_bytree'] = 0.8
    param['subsample'] = 0.8

    # gpu
    if args.gpu:
        param['gpu_id'] = int(args.gpu)
    #     param['tree_method'] = 'gpu_exact'
        param['tree_method'] = 'gpu_hist'
        param['max_bin'] = 1024
        param['predictor'] = 'cpu_predictor'

    print('Starting training...')
    print('xgboost params:')
    print(param)

    watchlist = [(d_train, 'train'), (d_test, 'eval')]

    num_round = 2000

    bst = xgb.train(param, d_train, num_round, watchlist, early_stopping_rounds=20)

    bst.save_model(model_file)

    scores = bst.get_score()
    ivar = 1
    for k in sorted(scores, key=scores.get, reverse=True):
        print("%2d. %24s: %s" % (ivar, k, str(scores[k])))
        ivar = ivar + 1

    # testing
    del bst
    bst = xgb.Booster()
    bst.load_model(model_file)

    pred = bst.predict(d_test)
    y_val = d_test.get_label()

    plotROC(pred, y_val)
#    plt.savefig(os.path.join(model_dir, 'roc.pdf'))
    plt.ylim(0.8, 1)
    plt.savefig(os.path.join(model_dir, 'roc_zoom.pdf'))

else:
    output_dir = os.path.dirname(args.output)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print('Start loading files...')
    test_filelist = glob.glob(args.input)
    d_test, df, y_truth = make_dmatrix(test_filelist, predict=True)
    print('Using testing files with %d events:\n%s\n...' % (d_test.num_row(), '\n'.join(test_filelist)))

    # setup parameters for xgboost
    bst = xgb.Booster()
    bst.load_model(model_file)

    preds = bst.predict(d_test)

    plotROC(preds, y_truth)
    # plt.savefig(os.path.join(model_dir, 'roc.pdf'))
    plt.ylim(0.8, 1)
    plt.savefig(os.path.join(output_dir, 'roc_zoom.pdf'))

    df['score_top'] = preds

    print('Write prediction file to %s' % args.output)
    import pandas as pd
    pandas_df = pd.DataFrame(df)
#     pandas_df.to_hdf(args.output.rsplit('.', 1)[0] + '.h5', 'Events', format='table')

    from root_numpy import array2root
    array2root(pandas_df.to_records(index=False), filename=args.output.rsplit(
        '.', 1)[0] + '.root', treename='Events', mode='RECREATE')
