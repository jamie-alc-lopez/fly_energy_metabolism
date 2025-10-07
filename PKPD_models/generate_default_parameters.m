%This script generates a structure of default parameters for the LG PKPD
%model. By default, it simulates the linear reaction network

clear;clc

%Reaction coefficients
par.kox_max = 10;
par.kaG = 100;
par.kaL = 100;
par.kgly = 1.5;
par.ldh_knockdown_frac = 1;
par.kldhf = 1;
par.kldhr = 0.01; 

%Gut dilution rate and central compartment clearance rates
par.delta = 0.5; 
par.keL = 0.1; 
par.keP = 0.1; 
par.keG = 0.1; 

%Diet parameters
par.Gdiet = 0.5;
par.Ldiet = 0.5;

%Parameters for OXPHOS regulation
par.rc = 0.5;
par.w = 0.13;
par.eps = 1e-7; 
par.central_ratio_reg = false; %Turn on or off reg
par.gut_ratio_reg = false;
par.kox_basal_frac = 0.01;

%Energy yield parameters
par.gly_yield = 2; %Per glucose glycolysis yield
par.ox_yield = 17; %Per pyruvate OXPHOS yield

%Energy balance model parameters, arbitrary
par.phi_O = 30;
par.max_phi_I = 19;
par.mu_E = 10;

save('default_parameters.mat',"par")