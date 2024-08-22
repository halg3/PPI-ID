# This script can take a directory as a part of a command-line argument to extract PAE, ipTM and pTM information.
# For example:
# > python3 alphafold_output_analysis.py /path/to/directory/containing/alphafold/output
# Some of this code was adopted from Lemal. (2022). "Explained: how to plot the prediction quality metrics with AlphaFold2," BioStrand, 
# https://blog.biostrand.ai/explained-how-to-plot-the-prediction-quality-metrics-with-alphafold2.

import os
import pickle
import json
import numpy as np
import matplotlib.pyplot as plt
import argparse

# Necessary functions
def get_pae_plddt(model_dicts):
    out = {}
    for i, d in enumerate(model_dicts):
        out[f'model_{i+1}'] = {'plddt': d['plddt'], 'pae': d['predicted_aligned_error']}
    return out

def generate_output_images(feature_dict):
    msa = feature_dict['msa']
    seqid = (np.array(msa[0] == msa).mean(-1))
    seqid_sort = seqid.argsort()
    non_gaps = (msa != 21).astype(float)
    non_gaps[non_gaps == 0] = np.nan
    final = non_gaps[seqid_sort] * seqid[seqid_sort, None]
    return final

# Class definition
class ARG:
    def __init__(self, repo):
        self.input_dir = repo
        self.output_dir = repo
        self.name = repo

def collective_pae(args):
    file_path = os.path.join(args.input_dir, 'features.pkl')
    feature_dict = pickle.load(open(file_path, 'rb'))

    model_dicts = [pickle.load(open(f'{args.input_dir}/result_model_{f}_multimer_v3_pred_0.pkl', 'rb'))
                   for f in range(1, 6)]

    pae_plddt_per_model = get_pae_plddt(model_dicts)
    generate_output_images(feature_dict)

    num_models = len(model_dicts)
    plt.figure(figsize=(3 * num_models, 2), dpi=100)
    for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
        plt.subplot(1, num_models, n + 1)
        plt.title(model_name)
        plt.imshow(value["pae"], label=model_name, cmap="bwr", vmin=0, vmax=30)
        plt.colorbar()
    plt.savefig(f"{args.output_dir}/collective_PAE.pdf")
    plt.close()

    return model_dicts

def best_pae(args):
    with open(os.path.join(args.input_dir, "ranking_debug.json"), 'r') as f:
        ranking_dict = json.load(f)

    best_model = ranking_dict["order"][0]

    file_path = os.path.join(args.input_dir, 'features.pkl')
    feature_dict = pickle.load(open(file_path, 'rb'))

    file_name = "result_" + best_model + ".pkl"
    file_path2 = os.path.join(args.input_dir, file_name)
    model_dict = [pickle.load(open(file_path2, 'rb'))]

    pae_plddt = get_pae_plddt(model_dict)
    generate_output_images(feature_dict)

    num_models = len(model_dict)
    plt.figure(figsize=(3 * num_models, 2), dpi=100)
    for n, (model_name, value) in enumerate(pae_plddt.items()):
        model_numb = best_model[0:7]
        prediction_numb = best_model[20:26]
        model_name = f'{model_numb}_{prediction_numb}'
        plt.title(model_name)
        plt.imshow(value["pae"], cmap="bwr", vmin=0, vmax=30)
        plt.colorbar()
    plt.savefig(f"{args.output_dir}/best_PAE.pdf")

    return best_model

def collective_iptm(args, model_dicts):
    iptm_scores = []
    for model in model_dicts:
        iptm_scores.append(float(model['iptm']))

    file_path_txt = os.path.join(args.output_dir, 'collective_iptm.txt')
    with open(file_path_txt, 'w') as file:
        file.write('\n'.join(map(str, iptm_scores)))

def best_iptm(args, best_model):
    best_model_dict = [pickle.load(open(f'{args.input_dir}/result_{best_model}.pkl', 'rb'))]

    for item in best_model_dict:
        if 'iptm' in item:
            iptm_value = [float(item['iptm'])]
            break

    file_path_txt = os.path.join(args.output_dir, 'best_iptm.txt')
    with open(file_path_txt, 'w') as file:
        file.write('\n'.join(map(str, iptm_value)))

def best_ptm(args, best_model):
    best_model_dict = [pickle.load(open(f'{args.input_dir}/result_{best_model}.pkl', 'rb'))]

    for item in best_model_dict:
        if 'ptm' in item:
            ptm_value = [float(item['ptm'])]
            break

    
    file_path_txt = os.path.join(args.output_dir, 'best_ptm.txt')
    with open(file_path_txt, 'w') as file:
        file.write('\n'.join(map(str, ptm_value)))

def main(repo):
    for r in repo:
        args = ARG(r)

        model_dicts = collective_pae(args)  # Generates PAE plots for all 5 models.
        best_model = best_pae(args)  # Generates PAE plot for best model
        collective_iptm(args, model_dicts)  # Fetches ipTM scores for all 5 models
        best_iptm(args, best_model)  # Fetches ipTM score for best model.
        best_ptm(args, best_model) # Fetches pTM score for best model.

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process directories containing AF output.")
    parser.add_argument("repo", nargs='+', help="Directories containing AF output")
    args = parser.parse_args()
    main(args.repo)
