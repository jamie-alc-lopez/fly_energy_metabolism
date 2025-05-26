#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 11:03:04 2024

@author: Jamie Alcira Lopez
"""

import pandas as pd
import operator
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns
from sklearn.decomposition import PCA

#Utility functions for analyzing metabolomic data
#Most of these functions are built to accomodate simultaneous processing of 
#negative and positive data

#This function imports metabolomics data, returning a manifest and two dictionaries
#The first dictionary contains the feature properties, while the second contains
#the peak areas in each sample
def import_metabolomics(manifest_file,data_file):
    #Load manifest with samples as index
    manifest = pd.read_excel(manifest_file,index_col=None)

    #Make unified sample name from negative sample name
    sample_list = [name.replace("_Neg_","_") for name in manifest.neg_name]
    manifest.insert(0,"sample_name",sample_list)
    manifest.set_index('sample_name',inplace=True)

    #Load all sheets in data, skipping the first two rows and using alignment ID as
    #index
    data_dict = pd.read_excel(data_file,index_col=1,skiprows=2,sheet_name=None)
    feature_dict = dict.fromkeys(data_dict.keys())
    peak_dict = dict.fromkeys(data_dict.keys())

    #Loop through all the metab data and get dictionaries of feature and peak dfs
    for key, ms_df in data_dict.items():
        
        #Switches between pos and neg names so right indexing is used,
        if ("Neg" in key) & ("Pos" not in key):
            ms_names = manifest.neg_name
        elif ("Neg" not in key) & ("Pos" in key):
            ms_names = manifest.pos_name
        else: 
            raise ValueError("Invalid dataset name; cannot determine run polarity.")
                        
        #Extract dataframe of feature information  
        final_feature_index = ms_df.columns.to_list().index("Adduct type")
        feature_info = ms_df.iloc[:, 0:(final_feature_index+1)] 

        #Extract feature peak AUC data corresponding to manifest names
        peak_data = ms_df.iloc[:,(final_feature_index+1):]
        peak_data = peak_data.T
        peak_data = peak_data.loc[ms_names]
        
        #Exchange polarity-specific index for universal sample names
        peak_data.insert(0,"sample_name",sample_list)
        peak_data.set_index('sample_name',inplace=True)
        
        #Remove features that are zero across all samples in the dataset
        nonzero_feature_IDs = peak_data.columns[peak_data.sum() > 0]
        peak_data = peak_data.loc[:,nonzero_feature_IDs]
        feature_info = feature_info.loc[nonzero_feature_IDs]
        
        #Add these to the dictionaries
        peak_dict[key] = peak_data
        feature_dict[key] = feature_info

    return manifest, feature_dict, peak_dict


#This function normalizes a dictionary of peak data by the specified
#normalization and pseudocount methods
def normalize_areas(feature_dict,peak_dict,norm_mode,pseudo_mode,**kwargs):
    
    #Initialize empty IS normalized dict
    peak_dict_n = dict.fromkeys(peak_dict.keys())

    #Loop through datasets
    for key in feature_dict.keys():
        
        feature_info = feature_dict[key]
        peak_data = peak_dict[key]
        
        #Add specified pseudocount
        if pseudo_mode == "min":
            pseudocount = peak_data[peak_data>0].min().min()
        elif pseudo_mode == "constant":
            pseudocount = kwargs.get('pseudocount',None)
        peak_data = peak_data + pseudocount
        
        #Identify the IS features
        IS_IDs = feature_info.index[feature_info.INCHIKEY == "Internal Standard"].to_list()
        
        if norm_mode == "mean_IS":
            norm_factor = peak_data.loc[:,IS_IDs].mean(axis = 1)
            peak_dict_n[key] = peak_data.div(norm_factor,axis=0)
        elif norm_mode == "none":
            peak_dict_n[key] = peak_data
            
    return peak_dict_n
        
#This function averages all peak abundances and returns a compressed manifest,
#and data on the averages and standard deviations of all peaks.
def average_replicates(manifest,peak_dict):
    
    #Get condition groupings, make new manifest
    g = manifest.groupby('condition')
    manifest_av = manifest.drop_duplicates(subset=['condition'])
    manifest_av.set_index('condition',inplace=True)
    manifest_av = manifest_av.drop(['neg_name','rep'],axis='columns')
    
    #Make empty new dictionaries
    peak_dict_av = dict.fromkeys(peak_dict.keys())
    peak_dict_std = dict.fromkeys(peak_dict.keys())
    
    #Loop through datasets
    for key, peak_data in peak_dict.items():
        peak_columns = peak_data.columns
        peak_index = manifest_av.index
        peak_dict_av[key] = pd.DataFrame(0,index=peak_index, columns=peak_columns)
        peak_dict_std[key] = pd.DataFrame(0,index=peak_index, columns=peak_columns)
        
        #Loop through conditions
        for condition in peak_dict_av[key].index:
            cond_samples = g.get_group(condition).index
            cond_peaks = peak_data.loc[cond_samples]
            peak_dict_av[key].loc[condition] = cond_peaks.mean(axis = 0)
            peak_dict_std[key].loc[condition] = cond_peaks.std(axis = 0)
    
    return manifest_av, peak_dict_av, peak_dict_std


#This function divides the peak AUCs to the average of a sample set
def norm_to_sample_set(peak_dict,sample_set):
    peak_dict_fc = dict.fromkeys(peak_dict.keys())
    for key, peak_data in peak_dict.items():
        if len(sample_set) > 1:
            norm_factor = peak_data.loc[sample_set].mean(axis = 0).values
        else:
            norm_factor = peak_data.loc[sample_set].values 
    
        peak_dict_fc[key] = peak_data.div(norm_factor,axis='columns')
    
    return peak_dict_fc


#This function returns a list of samples matching specific manifest properties
def subset_manifest(manifest,property_dict,comparison_dict):
    
    total_index = [True for i in range(len(manifest.index))]
    
    for key, value in property_dict.items():
        
        #If multiple inputs for a given property, take the "or" operation
        if isinstance(value, list):
            key_index = [False for i in range(len(manifest.index))]
            for subvalue in value:
                i_index = (manifest.loc[:,key] == subvalue).to_list()
                key_index = [a | b for a,b in zip(key_index,i_index)]
        else:    
            key_index = (manifest.loc[:,key] == value).to_list()
        
        #Negate if comparison is an exclusion
        if not comparison_dict[key]:
            key_index = list(map(operator.not_,key_index))
             
        #Update total index with this comparison
        total_index = [a & b for a,b in zip(total_index,key_index)]
    
    #Get final sample and manifest from comparison
    matching_samples = manifest.index[total_index].to_list()
    matching_manifest = manifest.loc[matching_samples]
    
    return matching_samples, matching_manifest
        

#This function takes in a set of sample names and returns the mean of the samples
def average_samples(sample_names,peak_df):
    sample_average = peak_df.loc[sample_names].mean()
    return sample_average
    
    
#This takes in a set of samples, an anchor sample, and a df of peak data and 
#generates a heatmap 
def generate_peak_heatmap(sample_names,anchor_sample,peak_df,label_names,title,
                          clabel='log10 normalized AUC',fontsize=12,cmap="BrBG",
                          cbarticks=None):
    
    peak_df = peak_df.loc[sample_names]
    peak_df.sort_values(axis=1,by=anchor_sample,ascending=False,inplace=True)
    plt_df = np.log10(peak_df)
    fig, ax = plt.subplots(figsize=(6.4,6.4))
    
    #Find symmetric points for vmin and vmax
    vmax = np.ceil(np.max((np.abs(plt_df.max().max()),np.abs(plt_df.min().min()))))
    vmin = -vmax
    ax = sns.heatmap(plt_df,cmap=cmap,
                     cbar_kws={'label': clabel,"ticks":cbarticks},
                     vmin = vmin, vmax = vmax)
    ax.set_yticks(np.array(range(len(label_names))) + 0.5)
    ax.set_xticks([])
    ax.set_yticklabels(label_names, fontsize=fontsize)
    ax.set_facecolor('#5ad9db')   

    plt.title(title,fontsize=fontsize)
    plt.xlabel('Metabolomic Features',fontsize=fontsize)
    plt.ylabel('')

    return fig, plt_df



#This mean-centers and standardizes variance of a df, returning a numpy array
def run_metab_pca(df):

    X = df.to_numpy()
    std_X = X.std(axis=0)
    Y = X[:,std_X > 0]
    std_Y = std_X[std_X>0]
    mean_Y = Y.mean(axis=0)
    Y = Y - mean_Y
    Y = Y/std_Y
    
    pca = PCA(n_components=4, svd_solver='full')
    pca.fit(Y)
    Y_transform = pca.transform(Y)

    return pca, Y_transform     


#Makes a PCA plot
def plot_pca(pca_manifest,pca,Y_transform,condition_colors,leg_labels,fontsize):
    for condition,color,alpha in condition_colors:
        
        cond_ind = pca_manifest['condition'] == condition
        plt.scatter(Y_transform[cond_ind,0],Y_transform[cond_ind,1],c=color,
                    alpha = alpha)
    plt.legend(leg_labels)

    plt.xlabel("PC1 ("+
               str(int(np.round(pca.explained_variance_ratio_[0],2)*100))+"%)",
               fontsize=fontsize)
    plt.ylabel("PC2 ("+
               str(int(np.round(pca.explained_variance_ratio_[1],2)*100))+"%)",
               fontsize=fontsize)
    plt.xticks([])
    plt.yticks([])

