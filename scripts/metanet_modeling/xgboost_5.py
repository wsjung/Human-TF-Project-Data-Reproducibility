import pandas as pd
import numpy as np
import os
import argparse
from sklearn.model_selection import cross_val_score
from hyperopt import fmin, tpe, hp, Trials, STATUS_OK
from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score, log_loss, average_precision_score, roc_curve, precision_recall_curve, auc, f1_score
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="XGBoost classifier train/test")
parser.add_argument("--tissue", required=True, help="tissue")
parser.add_argument("--base_path", required=True, help="base path")
parser.add_argument("--model", required=True, help="model (e.g. excl_perturbation)")
parser.add_argument("--binding_threshold", type=int, default=10, help="binding threshold in % (default=10)")
parser.add_argument("--n_cv_folds", type=int, default=10, help="number of CV folds")
parser.add_argument("--seed", type=int, default=42, help="random seed")
parser.add_argument("--output_subdir", type=str, default="xgboost", help="output subdirectory name")
args = parser.parse_args()

# set random seed
np.random.seed(seed=args.seed)

########################
### HYPEROPT CONFIGS ###
########################
# hyperopt xgboost parameter search space
space = {
    'max_depth': hp.quniform('max_depth', 1, 10, 1),
    'min_child_weight': hp.quniform('min_child_weight', 1, 5, 1),
    'subsample': hp.uniform('subsample', 0.5, 1.0),
    'colsample_bytree': hp.uniform('colsample_bytree', 0.5, 1.0),
    'gamma': hp.uniform('gamma', 0, 10),
    'reg_lambda': hp.uniform('reg_lambda', 0, 1),
    'learning_rate': hp.loguniform('learning_rate', np.log(0.01), np.log(0.3)),
    'n_estimators': hp.quniform('n_estimators', 20, 200, 20),
    'objective': 'binary:logistic',
    'eval_metric': 'logloss'
}


# paths
run_path = os.path.join(args.base_path, "results", args.tissue, f"{args.binding_threshold}.0_threshold_{args.model}")
input_cv_path = os.path.join(run_path, "CV_folds")
output_path = os.path.join(run_path, args.output_subdir)
output_cv_path = os.path.join(output_path, "CV")

# create dirs
if not os.path.exists(output_cv_path):
    os.makedirs(output_cv_path)

###############
### CV runs ###
###############
df_preds_list = []
for fold in range(args.n_cv_folds):
    print(f"fold: {fold}")

    # load CV data
    df_train = pd.read_csv(os.path.join(input_cv_path, f"fold{fold}_train_data.txt"), sep="\t")
    df_test = pd.read_csv(os.path.join(input_cv_path, f"fold{fold}_test_data.txt"), sep="\t")

    # split features vs response
    df_X_train = df_train.drop(columns=['TF','GENE','LABEL'])
    df_Y_train = df_train['LABEL']
    df_X_test = df_test.drop(columns=['TF','GENE','LABEL'])
    df_Y_test = df_test['LABEL']

    # class (im)balance
    class_balance_ratio = (len(df_Y_train) - sum(df_Y_train)) / sum(df_Y_train)
    print(f"class balance ratio: {class_balance_ratio}")
    # set scale_pos_weight
    space['scale_pos_weight'] = class_balance_ratio


    # hyperopt objective
    def objective(params):
        # Convert hyperparameters to integer where needed
        params['max_depth'] = int(params['max_depth'])
        params['min_child_weight'] = int(params['min_child_weight'])
        params['n_estimators'] = int(params['n_estimators'])

        # Initialize the model with the given hyperparameters
        model = XGBClassifier(**params)

        # cross-validation (use average_precision -- AUPRC deals better with class imbalances)
        cv_score = cross_val_score(model, df_X_train, df_Y_train, cv=5, scoring='average_precision', n_jobs=5).mean()

        # We want to minimize the negative accuracy (maximize accuracy)
        loss = -cv_score

        return {'loss': loss, 'status': STATUS_OK}

    # hyperopt tuning
    print("hyperparameter tuning with hyperopt")
    trials = Trials()
    best = fmin(fn = objective,
            space = space,
            algo = tpe.suggest,
            max_evals = 50,
            trials = trials)
    best['max_depth'] = int(best['max_depth'])
    best['min_child_weight'] = int(best['min_child_weight'])
    best['n_estimators'] = int(best['n_estimators'])
    df_best = pd.DataFrame(best.items(), columns=['param','value'])
    print("Best hyperparameters found: ", best)
    # write best params
    df_best.to_csv(os.path.join(output_cv_path, f"fold{fold}_hyperparams.tsv"), sep="\t", index=False)

    # train final model with best hyperparam values
    print("training final model")
    final_model = XGBClassifier(**best)
    final_model.fit(df_X_train, df_Y_train)

    # make predictions
    print("predicting")
    train_preds = final_model.predict(df_X_train)
    train_probs_class = final_model.predict_proba(df_X_train)
    train_probs = train_probs_class[:,1]

    test_preds = final_model.predict(df_X_test)
    test_probs_class = final_model.predict_proba(df_X_test)
    test_probs = test_probs_class[:,1]

    # evaluate
    train_accuracy = accuracy_score(df_Y_train, train_preds)
    train_logloss = log_loss(df_Y_train, train_probs)
    train_f1_score = f1_score(df_Y_train, train_preds)
    print(f'Training Accuracy: {train_accuracy}')
    print(f'Training Log Loss: {train_logloss}')
    print(f'Training F1 score: {train_f1_score}')

    test_accuracy = accuracy_score(df_Y_test, test_preds)
    test_logloss = log_loss(df_Y_test, test_probs)
    test_f1_score = f1_score(df_Y_test, test_preds)
    print(f'Test Accuracy: {test_accuracy}')
    print(f'Test Log Loss: {test_logloss}')
    print(f'Test F1 score: {test_f1_score}')

    # record predictions
    df_train_preds = df_train[['TF','GENE']].copy()
    df_train_preds['pred_prob'] = train_probs
    df_train_preds['pred'] = train_preds

    df_test_preds = df_test[['TF','GENE']].copy()
    df_test_preds['pred_prob'] = test_probs
    df_test_preds['pred'] = test_preds
    df_test_preds = df_test_preds.reset_index(drop=True)
    df_preds_list.append(df_test_preds)

    # plot auroc
    print("plotting AUROC")
    fpr_train, tpr_train, _ = roc_curve(df_Y_train, train_preds)
    roc_auc_train = auc(fpr_train, tpr_train)
    fpr_test, tpr_test, _ = roc_curve(df_Y_test, test_preds)
    roc_auc_test = auc(fpr_test, tpr_test)

    plt.figure(figsize=(10,10))
    plt.plot(fpr_train, tpr_train, color='darkorange', lw=2, label=f"ROC train (AUC={roc_auc_train:.2f})")
    plt.plot(fpr_test, tpr_test, color='blue', lw=2, label=f"ROC test (AUC={roc_auc_test:.2f})")
    plt.plot([0,1], [0,1], lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f"{args.tissue}\nFold {fold} AUROC")
    plt.legend()
    plt.savefig(os.path.join(output_cv_path, f"fold{fold}_auroc.png"), format="png", bbox_inches="tight")

    # plot auprc
    print("plotting AUPRC")
    precision_train, recall_train, _ = precision_recall_curve(df_Y_train, train_probs)
    prc_auc_train = auc(recall_train, precision_train)
    precision_test, recall_test, _ = precision_recall_curve(df_Y_test, test_probs)
    prc_auc_test = auc(recall_test, precision_test)

    plt.figure(figsize=(10,10))
    plt.plot(recall_train, precision_train, color='darkorange', lw=2, label=f"PRC train (AUC={prc_auc_train:.2f})")
    plt.plot(recall_test, precision_test, color='blue', lw=2, label=f"PRC test (AUC={prc_auc_test:.2f})")
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(f"{args.tissue}\nFold {fold} AUPRC")
    plt.legend()
    plt.savefig(os.path.join(output_cv_path, f"fold{fold}_auprc.png"), format="png", bbox_inches="tight")

    print("saving")
    # write CV preds
    df_test_preds = df_test_preds.sort_values(by='pred_prob', ascending=False)
    df_test_preds.to_csv(os.path.join(output_cv_path, f"fold{fold}_test_preds.tsv"), sep="\t", index=False)
    # save object for downstream analysis (e.g. SHAP)
    final_model.save_model(os.path.join(output_cv_path, f"fold{fold}_model.json"))

# write final preds
df_preds = pd.concat(df_preds_list)
df_preds = df_preds.sort_values(by='pred_prob', ascending=False)
df_preds.to_csv(os.path.join(output_path, f"{args.tissue}_xgboost_full.tsv"), sep="\t", index=False)
# tf-gene-pred_proba for evaluation
df_preds = df_preds.drop(columns=['pred'])
df_preds.to_csv(os.path.join(output_path, f"{args.tissue}_xgboost.tsv"), sep="\t", index=False, header=False)



print("DONE")
