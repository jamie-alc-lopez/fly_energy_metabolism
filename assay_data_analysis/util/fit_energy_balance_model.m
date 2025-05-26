function final_fit = fit_energy_balance_model(manifest)

%This function fits data from fly starvation experiments involving flies
%fed dilutions of a single type of media. The fit is performed in two
%stages, with first the mean being fit followed by the variance

% Set up regression on the mean survival time
ft_mu = fittype('model_mu_td(phi_I_base,phi_O,mu_E,f)',...
    'coefficients',{'phi_I_base','phi_O','mu_E'},'independent','f');
opts = fitoptions( 'Method', 'NonlinearLeastSquares');
opts.Lower = [0,0,1];
opts.StartPoint = [0.175 0.4,1];
opts.Upper = [1e5 1e5,1];

%Perform the fit and get confidence intervals
[fitresult_mu, gof_mu] = fit( manifest.nutrition, manifest.mu_td, ft_mu, opts);
c_ints_mu = confint(fitresult_mu); 
c_width_mu = [fitresult_mu.phi_I_base; fitresult_mu.phi_O] - c_ints_mu(1,1:2)';

%Extract results
final_fit.phi_I_base = fitresult_mu.phi_I_base;
final_fit.phi_O = fitresult_mu.phi_O;
final_fit.mu_E = fitresult_mu.mu_E;
final_fit.confint_phi_I_base = c_width_mu(1);
final_fit.confint_phi_O = c_width_mu(2);
final_fit.ratio_I_O = fitresult_mu.phi_I_base/fitresult_mu.phi_O; 
final_fit.mu_adj_rsquared = gof_mu.adjrsquare;

%Perform regression on the survival time variances
ft_sigma = fittype('model_sigma_td(sigma_E,phi_I_base,phi_O,f)',...
    'coefficients',{'sigma_E'},'independent',{'f'},'problem',{'phi_I_base','phi_O'});
opts = fitoptions( 'Method', 'NonlinearLeastSquares');
opts.Lower = 0;
opts.StartPoint = 0.1;
opts.Upper = 1e3;

%Get confidence intervals on variance params
[fitresult_sigma,gof_sigma] = fit( manifest.nutrition, manifest.sigma_td, ft_sigma, opts,...
    'problem',{fitresult_mu.phi_I_base,fitresult_mu.phi_O});
c_ints_sigma = confint(fitresult_sigma); 
c_width_sigma = [fitresult_sigma.sigma_E] - c_ints_sigma(1,:)';

final_fit.sigma_E = fitresult_sigma.sigma_E;
final_fit.confint_sigma_E = c_width_sigma;
final_fit.sigma_adj_rsquared = gof_sigma.adjrsquare;



end