function sigma_td = model_sigma_td(sigma_E,phi_I_base,phi_O,f)

%This function is the energy balance model prediction for the survival
%time standard deviation. It is used in fitting the energy balance model to
%assay data

sigma_td = sigma_E./(phi_O - f.*phi_I_base);

end