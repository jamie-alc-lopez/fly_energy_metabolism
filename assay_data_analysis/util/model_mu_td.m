function mu_td = model_mu_td(phi_I_base,phi_O,mu_E,f)

%This function is the energy balance model prediction for the mean survival
%time. It is used in fitting the energy balance model to assay data

mu_td = mu_E./(phi_O - f.*phi_I_base);

end