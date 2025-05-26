%% This script takes in data from the glc lac mixture 1x 2x experiment
% and estimates the parameters of an energy balance model

%We estimate the parameters for glucose, and then use these estimate the
%changes needed in energy intake, output, and storage needed to reproduce
%the effect of the other mixtures

clear;clc

%Load glc lac mixture moments
data = readtable('estimated_assay_moments/glc_lac_mixture_1x_2x_moments.xlsx');

%Specify media types
media_type_cell = {'glc','lac','mix'};

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
writetable(fit_table,'fit_results/glc_lac_mixture_1x_2x_fit.xlsx',...
    'WriteRowNames',true);


%% Estimate parameter changes needed for survival time change between glucose 
% vs. pure lactic acid and lactic acid/glucose at different caloric densities

%List pairs of conditions to analyze
cond_pair_cell = {{'female_62_5_glc','female_125_lac'},...
    {'female_62_5_glc','female_31_25_glc_62_5_lac'},...
    {'female_125_glc','female_250_lac'},...
    {'female_125_glc','female_62_5_glc_125_lac'}};

%Generate table to store energy delta estimates
variable_name_cell = {'cond1','cond2','t2_div_t1','fc_mu_E',...
    'fc_phi_I','fc_phi_O'};
energy_delta_table = cell(length(cond_pair_cell),length(variable_name_cell));
energy_delta_table = array2table(energy_delta_table,...
    'VariableNames',variable_name_cell);

%Extract index corresponding to pure glucose fit
glc_fit_ind = strcmp(fit_table.Properties.RowNames,'glc');

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
    phi_O = fit_table.phi_O(glc_fit_ind);
    if strcmp(energy_delta_table.cond1{i},'female_62_5_glc')
        phi_I = fit_table.phi_I_base(glc_fit_ind);
    elseif strcmp(energy_delta_table.cond1{i},'female_125_glc')
        phi_I = 2*fit_table.phi_I_base(glc_fit_ind);
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
writetable(energy_delta_table,'fit_results/glc_lac_mixture_1x_2x_parameter_changes.xlsx');