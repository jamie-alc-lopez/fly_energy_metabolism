function [dy,kox_eff,kgly_eff,energy_flux] = LG_PKPD_RHS(t,y)

%This function computes the RHS of the lactic acid and glucose PKPD model
%model with ratio-dependent OXPHOS regulation. 

global par

%Extract state variables
Ggut = y(1);
Lgut = y(2);
G = y(3);
L = y(4);
P = y(5);

%Compute effective OXPHOS activity
if par.central_ratio_reg
    glu_lac_ratio = (G+par.eps)/(L+par.eps);
    kox_eff =  par.kox_basal_frac*par.kox_max+(1-par.kox_basal_frac)...
        *par.kox_max.*exp(-(glu_lac_ratio - par.rc).^2./(2*par.w^2));
elseif par.gut_ratio_reg
    glu_lac_ratio = (Ggut+par.eps)/(Lgut+par.eps);
    kox_eff =  par.kox_basal_frac*par.kox_max+(1-par.kox_basal_frac)...
        *par.kox_max.*exp(-(glu_lac_ratio - par.rc).^2./(2*par.w^2));
else
    kox_eff = par.kox_basal_frac*par.kox_max;
end

%Compute effective LDH rates with knockdown
kldhf_eff = par.ldh_knockdown_frac*par.kldhf;
kldhr_eff = par.ldh_knockdown_frac*par.kldhr;

%Define effective glycolysis rate (placeholder for models extensions with
%glycolysis regulation
kgly_eff = par.kgly; 

%Define the reaction function (placeholder for model extensions with
%Hill reaction kinetics
rxn_fun = @(x,K) x;

%Compute RHS for non-redox compounds
dy = zeros(size(y));

%Gut glucose
dy(1) = par.delta*par.Gdiet - par.kaG*rxn_fun(Ggut) - par.delta*Ggut;

%Gut lactic acid
dy(2) = par.delta*par.Ldiet - par.kaL*rxn_fun(Lgut) - par.delta*Lgut;

%Central compartment glucose
dy(3) = par.kaG*rxn_fun(Ggut) - kgly_eff*rxn_fun(G) ...
    - par.CLG*G;

%Central compartment lactic acid
dy(4) = par.kaL*rxn_fun(Lgut) - kldhf_eff*rxn_fun(L) ...
    + kldhr_eff*rxn_fun(P) - par.CLL*L;

%Central compartment pyruvate
dy(5) = 2*kgly_eff*rxn_fun(G) - ...
    kox_eff*rxn_fun(P) + kldhf_eff*rxn_fun(L)...
    - kldhr_eff*rxn_fun(P) - par.CLP*P;

%Compute the energy flux
energy_flux = par.gly_yield*kgly_eff*rxn_fun(G) + par.ox_yield*kox_eff*rxn_fun(P);


end