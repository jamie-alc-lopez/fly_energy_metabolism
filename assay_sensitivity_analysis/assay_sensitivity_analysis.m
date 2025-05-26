%This script performs a sensitivity analysis of the starvation assay as a
%function of the dietary energy content. The energy balance model is
%parameterized with values from the "orig_spent_media" CDM experiment. We
%use a finite lifespan model 

clear;clc

%Parameters from fit
mu_E = 1;
phi_O = 0.31;

%Set maximimum lifespan, compute corresponding maximal intake energy flux
mu_td_max = 30; 
phi_I_max = (phi_O*mu_td_max - mu_E)/mu_td_max;

%Vector of input energy levels
n = 200;
nutrition_vec = linspace(0,2,n);
phi_I_vec = nutrition_vec.*phi_I_max;
phi_I_vec(phi_I_vec>phi_I_max) = phi_I_max;

%Compute survival times
mu_td = mu_E./(phi_O-phi_I_vec);

%Compute sensitivities
h = max(nutrition_vec)/(n+1);
assay_sensitivity = diff(mu_td)/h;
delta_assay_sensitivity = diff(assay_sensitivity)/h;
assay_sensitivity = [assay_sensitivity,NaN];
delta_assay_sensitivity = [delta_assay_sensitivity,NaN,NaN];

%Make table of sensitivity analysis
sensitivity_table = array2table([nutrition_vec',mu_td',...
    assay_sensitivity',delta_assay_sensitivity'],...
    "VariableNames",{'diet_energy','mu_td','assay_sensitivity',...
    'delta_assay_sensitivity'});

%Save table
writetable(sensitivity_table,'assay_sensitivity_results.xlsx',...
    'WriteVariableNames',true)
