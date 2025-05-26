function [mu,sigma,CV] = compute_td_moments(t,S,only_real)

%This script takes in a time vector (t) and survival vector from one vial
% (S) and computes the mean, standard deviation, and CV. Numerical
% integration via Simpson's method is used to compute integrals. 

%Truncate the data to eliminate NaNs
nan_ind = find(isnan(S))-1;
trunc_point = min([min(nan_ind),length(S)]);
S_eff = S(1:trunc_point);
t_eff = t(1:trunc_point);

%Normalize data to be survival fraction
if max(S_eff) > 1
    S_eff = S_eff/max(S_eff);
end

%Only compute moments if the vial timecourse is complete (i.e. all flies
%eventually died). 
if S_eff(end) == 0

    %Numerically integrate to get moments
    mu = simps(t_eff,S_eff);
    sigma = sqrt(2*simps(t_eff,t_eff.*S_eff) - mu^2);

    %Eliminate imaginary sigma estimates if "only_real" selected
    %These originate from marginally negative variance estimates
    if only_real && ~isreal(sigma)
        sigma = 0;
    end

    %Compute CV
    CV = sigma/mu;

%If timecoure is incomplete, return NaN
else
    mu = nan;
    sigma = nan;
    CV = nan;
end

end