import pandas as pd
import numpy as np
import os
import argparse
import sys


def eval_binding(df_network, df_label, num_tfs, max_avg_targets, outputfile):


    label_keys = set(df_label['key'])

    df_network['is_supported'] = df_network['key'].isin(label_keys).astype(int)
    df_network['cum_supported'] = df_network['is_supported'].cumsum()

    rows = []
    for avg in range(1,max_avg_targets+1):
        K = avg * num_tfs # number of edges
        supported = int(df_network.loc[K-1, "cum_supported"])
        perc_support = supported / K * 100
        rows.append({
            "avg_targets_per_tf": avg,
            "binding_support": perc_support
        })

    df_out = pd.DataFrame(rows)
    df_out.to_csv(outputfile, sep="\t", index=False, header=False)

def eval_binding_full(df_network, df_label, num_tfs, outputfile):

    label_keys = set(df_label['key'])
    df_network['is_supported'] = df_network['key'].isin(label_keys).astype(int)
    df_network['cum_supported'] = df_network['is_supported'].cumsum()

    perc_support = df_network.iloc[-1]['cum_supported'] / len(df_network) * 100

    df_out = pd.DataFrame({
        'threshold':['full'],
        'perc':[perc_support]
    })
    df_out.to_csv(outputfile, sep="\t", index=False, header=False)



def prep_args(network, label):

    df_network = pd.read_csv(network, sep="\t", header=None, names=['REGULATOR','TARGET','SCORE'])
    df_network = df_network.sort_values(by='SCORE', ascending=False).reset_index(drop=True)

    df_label = pd.read_csv(label, sep="\t")

    df_label = df_label[df_label['REGULATOR'].isin(df_network['REGULATOR'])]

    ntfs_network = df_network['REGULATOR'].nunique()
    ntfs_label = df_label['REGULATOR'].nunique()

    # assign keys
    df_network['key'] = df_network['REGULATOR'] + '-' + df_network['TARGET']
    df_label['key'] = df_label['REGULATOR'] + '-' + df_label['TARGET']

    return df_network, df_label, ntfs_network


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="evaluates binding metric")
    parser.add_argument("--network", required=True, help="input network")
    parser.add_argument("--label", required=True, help="binding label")
    parser.add_argument("--max_avg_targets", type=int, default=-1, help="max average targets per TF threshold")
    parser.add_argument("--output", required=True, help="output filepath")
    parser.add_argument("--full_network", action='store_true', help="use to evaluate the full network, instead of at thresholds")
    args = parser.parse_args()

    if not args.full_network and args.max_avg_targets == -1:
        print("ERROR: argument --max_avg_targets <int> not specified without --full_network argument.")
        sys.exit(1)

    df_network, df_label, num_tfs = prep_args(args.network, args.label)

    if args.full_network:
        print("> evaluating full network")
        eval_binding_full(df_network, df_label, num_tfs, args.output)
    else:
        print("> evaluating network at thresholds")
        eval_binding(df_network, df_label, num_tfs, args.max_avg_targets, args.output)
    print("DONE")
    sys.exit(0)


