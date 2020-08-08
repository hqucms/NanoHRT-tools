import numpy as np
import xgboost as xgb


class XGBHelper:

    def __init__(self, model_file, var_list):
        self.bst = xgb.Booster(params={'nthread': 1}, model_file=model_file)
        self.var_list = var_list
        print('Load XGBoost model %s, input variables:\n  %s' % (model_file, str(var_list)))

    def eval(self, inputs):
        dmat = xgb.DMatrix(np.array([[inputs[k] for k in self.var_list]]), feature_names=self.var_list)
        return self.bst.predict(dmat)[0]


class XGBEnsemble:

    def __init__(self, model_files, var_list):
        self.bst_list = [xgb.Booster(params={'nthread': 1}, model_file=f) for f in model_files]
        self.var_list = var_list
        print('Load XGBoost models:\n  %s, \ninput variables:\n  %s' % ('\n  '.join(model_files), str(var_list)))

    def eval(self, inputs, model_idx=None):
        dmat = xgb.DMatrix(np.array([[inputs[k] for k in self.var_list]]), feature_names=self.var_list)
        if model_idx is not None:
            return self.bst_list[model_idx].predict(dmat)[0]
        else:
            preds = [bst.predict(dmat)[0] for bst in self.bst_list]
            return sum(preds) / len(self.bst_list)
