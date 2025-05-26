clear;clc

%This script estimates the moments (mean and variance) of all assay data
%within the "raw_data" folder and saves the resulting tables in the
%"estimated_assay_moments" folder.

%Get list of all raw data files
data_folder = 'raw_assay_data';
data_files = dir([data_folder,'/*.xlsx']);

%Specify output folder
output_folder = 'estimated_assay_moments';

%In rare edge cases the numerical estimator leads to imaginary
%moment estimates. Setting "only_real" to "true" removes the imaginary
%component in these cases
only_real = true;

%For fun, track total number of flies used in assays
total_flies = 0; 

%Loop through datasets
for i = 1:length(data_files)

    %Get dataset name
    file_name = data_files(i).name;
    dataset_name = strrep(file_name,'_data.xlsx','');

    %Avoid running on temporary files
    if ~contains(dataset_name,'~$')

        %Load raw data and associated manifest
        data_i = readtable([data_folder,'/',file_name],'Sheet','survival_data');
        manifest_i = readtable([data_folder,'/',file_name],'Sheet','manifest');

        %Extract time vector
        t_i = data_i.time;

        %Loop through the different vials
        for j = 1:size(manifest_i,1)

            vial_name = matlab.lang.makeValidName(manifest_i.name{j});

            %Compute and store moments
            S_ij = data_i{:,vial_name};
            [manifest_i.mu_td(j), manifest_i.sigma_td(j),manifest_i.CV(j)] ...
                = compute_td_moments(t_i,S_ij,only_real);

        end

        %Export the estimated moments
        writetable(manifest_i,[output_folder,'/',dataset_name,'_moments.xlsx'],...
            'WriteVariableNames',true,'WriteRowNames',true);

        total_flies = total_flies + size(manifest_i,1)*10;

    end

end


%% Perform post-switch death analysis for priming experiment

%Get list of all raw data files
data_folder = 'raw_assay_data';

%Specify output folder
output_folder = 'estimated_assay_moments';

%Get dataset name
file_name = 'priming_data.xlsx';
dataset_name = strrep(file_name,'_data.xlsx','');

%Load raw data and associated manifest
data = readtable([data_folder,'/',file_name],'Sheet','survival_data');
manifest = readtable([data_folder,'/',file_name],'Sheet','manifest');

%Remove starvation condition as none survive at switch time
data = data(:,~contains(data.Properties.VariableNames,'stv'));
manifest = manifest(~strcmp(manifest.condition,'stv'),:);

%Extract time vector
t = data.time;
switch_time = 5; 
switch_ind = find(t == switch_time);
t_eff = t(switch_ind:end);
t_eff = t_eff - switch_time;

%In rare edge cases the numerical estimator leads to imaginary
%moment estimates. Setting "only_real" to "true" removes the imaginary
%component in these cases
only_real = true;

%Loop through the different vials
for j = 1:size(manifest,1)

    vial_name = matlab.lang.makeValidName(manifest.name{j});

    %Compute and store moments
    S_j = data{:,vial_name};
    S_j_eff = S_j(switch_ind:end);
    [manifest.mu_td(j), manifest.sigma_td(j),manifest.CV(j)] ...
        = compute_td_moments(t_eff,S_j_eff,only_real);

end

%Export the estimated moments
writetable(manifest,[output_folder,'/',dataset_name,'_post_switch_moments.xlsx'],...
    'WriteVariableNames',true,'WriteRowNames',true);


