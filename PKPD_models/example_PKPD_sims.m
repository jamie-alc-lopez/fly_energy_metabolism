%This script produces the simulation results for the linear,
%OXPHOS-regulation and LDH knockdown model figures

%% Generate linear model simulation data

clear;clc

global par

load("default_parameters.mat")

%Initialize state variable storage
y0 = ones(5,1);
tspan = [0,1e10];
influx_total_vec = [1 2];
n = 200;
frac_glu = linspace(0,1,n);

for i = 1:length(influx_total_vec)
    for j = 1:n
        par.Gdiet = influx_total_vec(i)*frac_glu(j);
        par.Ldiet = influx_total_vec(i)*2*(1-frac_glu(j));
        [yend{i,j},energy_flux(i,j),kox_eff(i,j),~,~,rel_dy{i,j}] ...
            = run_LG_PKPD(y0,tspan);
    end
end

%Energy balance model
lifespan = par.mu_E./(par.phi_O-energy_flux);
lifespan(energy_flux > par.max_phi_I) = par.mu_E./(par.phi_O-par.max_phi_I);

%Generate lifespan table for figure
lifespan_table = ...
    array2table([frac_glu',lifespan',energy_flux(1,:)',energy_flux(2,:)'],...
    'VariableNames',{'Fraction glucose','1x lifespan','2x lifespan',...
    '1x energy flux','2x energy flux'});

%Save file
writetable(lifespan_table,...
    'sim_results/linear_LG_PKPD_lifespan.xlsx',...
    'WriteVariableNames',true)


%% Generate OXPHOS ratio regulation simulation data

clear;clc

global par

load("default_parameters.mat")

%Activate ratio regulation based on the gut
par.gut_ratio_reg = 1;

%Initialize state variable storage
y0 = ones(5,1);
tspan = [0,1e10];
influx_total_vec = [1 2];
n = 200;
frac_glu = linspace(0,1,n);

for i = 1:length(influx_total_vec)
    for j = 1:n
        par.Gdiet = influx_total_vec(i)*frac_glu(j);
        par.Ldiet = influx_total_vec(i)*2*(1-frac_glu(j));
        [yend{i,j},energy_flux(i,j),kox_eff(i,j),~,~,rel_dy{i,j}] ...
            = run_LG_PKPD(y0,tspan);
    end
end

%Energy balance model
lifespan = par.mu_E./(par.phi_O-energy_flux);
lifespan(energy_flux > par.max_phi_I) = par.mu_E./(par.phi_O-par.max_phi_I);

%Generate lifespan table for figure
lifespan_table = ...
    array2table([frac_glu',lifespan'],'VariableNames',...
    {'Fraction glucose','1x lifespan','2x lifespan'});

writetable(lifespan_table,...
    'sim_results/ratio_regulation_LG_PKPD_lifespan.xlsx',...
    'WriteVariableNames',true)


%% Generate LDH knockdown simulation data

clear;clc

global par

load("default_parameters.mat")

%Activate ratio regulation based on the gut
par.gut_ratio_reg = 1;

%Initialize state variable storage
y0 = ones(5,1);
tspan = [0,1e10];
influx_total = 1;
ldh_kd_frac_vec = [1 0.5 0.1 0.01];
n = 200;
frac_glu = linspace(0,1,n);

for i = 1:length(ldh_kd_frac_vec)
    par.ldh_knockdown_frac = ldh_kd_frac_vec(i);
    for j = 1:n
        par.Gdiet = influx_total*frac_glu(j);
        par.Ldiet = influx_total*2*(1-frac_glu(j));
        [yend{i,j},energy_flux(i,j),kox_eff(i,j),~,~,rel_dy{i,j}] ...
            = run_LG_PKPD(y0,tspan);
    end
end

%Energy balance model
lifespan = par.mu_E./(par.phi_O-energy_flux);
lifespan(energy_flux > par.max_phi_I) = par.mu_E./(par.phi_O-par.max_phi_I);

%Generate lifespan table for figure
lifespan_table = ...
    array2table([frac_glu',lifespan'],'VariableNames',...
    {'Fraction glucose','0% KD','50% KD','90% KD','99% KD'});

writetable(lifespan_table,...
    'sim_results/LDH_knockdown_LG_PKPD_lifespan.xlsx',...
    'WriteVariableNames',true)
