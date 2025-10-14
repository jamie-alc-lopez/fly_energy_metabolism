#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 11:26:49 2024

@author: Jamie Alcira Lopez
"""

#This script performs the metabolomics analysis for the main text figures

import metabolomics_util as mu
import pandas as pd
import matplotlib.pyplot as plt
#import importlib

fontsize = 16

#%% Import data and perform initial processing
#importlib.reload(mu)

#Import data
manifest_file = "metabolomics_data/spent_media_metabolomics_manifest.xlsx"
data_file = "metabolomics_data/spent_media_metabolomics_data.xlsx"
manifest, feature_dict, peak_dict = mu.import_metabolomics(manifest_file,data_file)

#Normalize and pseudocount the samples
norm_mode = "none"
pseudo_mode = "min"
peak_dict_n = mu.normalize_areas(feature_dict,peak_dict,norm_mode,pseudo_mode)

#Average replicates
manifest_av, peak_dict_av, peak_dict_std = mu.average_replicates(manifest,peak_dict_n)

#Compute version of data that is only CDM fold changes
property_dict = {"media":'CDM',"microbe":'none'}
comparison_dict = {"media":True,"microbe":True}
CDM_samples, __ = mu.subset_manifest(manifest_av,property_dict,comparison_dict)
peak_dict_fc = mu.norm_to_sample_set(peak_dict_av,CDM_samples)

#%% PCA visualization of the three main spent media and CDM

#Get desired subset of samples
microbes = ["LpWF","Ap","LpWFAp","none"]
property_dict = {"media":'CDM',"microbe":microbes,"time":[0,48]}
comparison_dict = {"media":True,"microbe":True,"time":True}
pca_samples, pca_manifest = mu.subset_manifest(manifest,property_dict,comparison_dict)
pca_df = peak_dict_n['Neg'].loc[pca_samples].copy()

#Perform PCA, get components
pca, Y_transform = mu.run_metab_pca(pca_df)

#Plot the PCA to check
pca_fig = plt.figure(1, figsize=(3, 3))
condition_colors = [("CDM","#D1735A",1),("48hr_sCDM_LpWF","#8DC3D8",1),
          ("48hr_sCDM_Ap","#6AA99B",1),("48hr_sCDM_LpWFAp","#586992",1)] #CDM, Lp, Ap, LpAp
leg_labels = {"CDM","Lp","Ap","LpAp"}
mu.plot_pca(pca_manifest,pca,Y_transform,condition_colors,leg_labels,fontsize)

#Save the data
pca_df = pd.DataFrame(Y_transform,index=pca_samples,columns=["PC1","PC2","PC3","PC4"])
pca_df.to_excel("processed_data/spent_media_pca_data.xlsx")


#%% Data for heatmaps of fold-change vs. CDM 

#Set up data
microbes = ["LpWF","Ap","LpWFAp"]
dataset = 'Neg'
peak_df = peak_dict_fc[dataset]
property_dict = {"media":'CDM',"microbe":microbes,"time":48}
comparison_dict = {"media":True,"microbe":True,"time":True}
__, sample_manifest = mu.subset_manifest(manifest_av,property_dict,comparison_dict)
sample_manifest.sort_values(by = ['time'],inplace=True)
sample_names = ["48hr_sCDM_LpWF","48hr_sCDM_Ap","48hr_sCDM_LpWFAp"]
label_names = ["Lp","Ap","LpAp"]
anchor_sample = sample_names[0]

#Plot the heatmap to check
title = None
clabel='log10(fold change vs. CDM)'
cbarticks=(-8,-4,0,4,8)
heatmap_fig, plt_df = mu.generate_peak_heatmap(sample_names,anchor_sample,peak_df,
                               label_names,title,clabel,cmap='BrBG',
                               fontsize=fontsize)
heatmap_fig.set_figheight(2.4)

#Save the data
plt_df.to_excel('processed_data/CDM_fold_change_heatmap_data.xlsx')   

