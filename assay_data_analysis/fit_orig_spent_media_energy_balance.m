%% This script takes in data from the original spent media experiment (i.e.,
%CDM vs. spent media vs. starvation in female OrR flies) and estimates the
%parameters of an energy balance model

%We estimate the parameters for CDM, and then use these estimate the
%changes needed in energy intake, output, and storage needed to reproduce
%the effect of the spent media

%We use this data to estimate the CDM parameters, rather than the data from
%the multidilution experiment, in order to ensure maximal comparability 
%with the spent media data

clear;clc

%Load spent media moments
data = readtable('estimated_assay_moments/orig_spent_media_moments.xlsx');

%Specify media types
media_type_cell = {'cdm','lp','ap','lpap'};

%Generate table to store fit results
variable_name_cell = {'phi_I_base','phi_O','mu_E','sigma_E',...
    'confint_phi_I_base','confint_phi_O','confint_sigma_E','mu_adj_rsquared',...
    'sigma_adj_rsquared'};
fit_table = nan(length(media_type_cell),length(variable_name_cell));
fit_table = array2table(fit_table,'RowNames',media_type_cell,...
    'VariableNames',variable_name_cell);

%Loop through media conditions and fit models
for i = 1:length(media_type_cell)

    %Get media type data
    media_ind = strcmp(data.exp,media_type_cell{i}) | strcmp(data.exp,'stv');
    subdata = data(media_ind,:);

    %Perform fit
    efit = fit_energy_balance_model(subdata);

    %Store results
    fit_table.phi_I_base(i) = efit.phi_I_base;
    fit_table.phi_O(i) = efit.phi_O;
    fit_table.mu_E(i) = efit.mu_E;
    fit_table.sigma_E(i) = efit.sigma_E;
    fit_table.confint_phi_I_base(i) = efit.confint_phi_I_base;
    fit_table.confint_phi_O(i) = efit.confint_phi_O;
    fit_table.confint_sigma_E(i) = efit.confint_sigma_E;
    fit_table.mu_adj_rsquared(i) = efit.mu_adj_rsquared;
    fit_table.sigma_adj_rsquared(i) = efit.sigma_adj_rsquared;

end

%Export table
writetable(fit_table,'fit_results/orig_spent_media_fit.xlsx',...
    'WriteRowNames',true);


%% Estimate parameter changes needed for survival time change between cdm 
% and spent media at different caloric densities

%List pairs of conditions to analyze
cond_pair_cell = {{'cdm_1x','lpwf_scdm_1x'},...
    {'cdm_1x','ap_scdm_1x'},...
    {'cdm_1x','lpwf_ap_scdm_1x'},...
    {'cdm_2x','lpwf_scdm_2x'},...
    {'cdm_2x','ap_scdm_2x'},...
    {'cdm_2x','lpwf_ap_scdm_2x'}};

%Generate table to store energy delta estimates
variable_name_cell = {'cond1','cond2','t2_div_t1','fc_mu_E',...
    'fc_phi_I','fc_phi_O'};
energy_delta_table = cell(length(cond_pair_cell),length(variable_name_cell));
energy_delta_table = array2table(energy_delta_table,...
    'VariableNames',variable_name_cell);

%Extract index corresponding to CDM fit
cdm_fit_ind = strcmp(fit_table.Properties.RowNames,'cdm');

%Loop through the different condition pairs
for i = 1:length(cond_pair_cell)

    %Store condition names
    energy_delta_table.cond1{i} = cond_pair_cell{i}{1};
    energy_delta_table.cond2{i} = cond_pair_cell{i}{2};

    %Extract mean survival times in the two conditions
    t1 = mean(data.mu_td(strcmp(data.condition,energy_delta_table.cond1{i})));
    t2 = mean(data.mu_td(strcmp(data.condition,energy_delta_table.cond2{i})));

    %Get survival time fold changes
    energy_delta_table.t2_div_t1{i} = t2/t1;

    %Extract fitting parameters
    phi_O = fit_table.phi_O(cdm_fit_ind);
    if strcmp(energy_delta_table.cond1{i},'cdm_1x')
        phi_I = fit_table.phi_I_base(cdm_fit_ind);
    elseif strcmp(energy_delta_table.cond1{i},'cdm_2x')
        phi_I = 2*fit_table.phi_I_base(cdm_fit_ind);
    else
        error('Invalid base condition specified.')
    end
    phi_I_div_phi_O = phi_I/phi_O;

    %Compute changes required in energy metabolism parameters
    energy_delta_table.fc_mu_E{i} = energy_delta_table.t2_div_t1{i};

    energy_delta_table.fc_phi_O{i} = ...
        1/energy_delta_table.t2_div_t1{i} - ...
        phi_I_div_phi_O*(1/energy_delta_table.t2_div_t1{i}) + phi_I_div_phi_O;

    energy_delta_table.fc_phi_I{i} = ...
        -(1/phi_I_div_phi_O)*(1/energy_delta_table.t2_div_t1{i}) ...
        + 1/energy_delta_table.t2_div_t1{i} + (1/phi_I_div_phi_O);

end

%Export table
writetable(energy_delta_table,'fit_results/orig_spent_media_parameter_changes.xlsx');