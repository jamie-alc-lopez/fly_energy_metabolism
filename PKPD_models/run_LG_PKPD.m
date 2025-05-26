function [yend,energy_flux,kox_eff,tout,yout,rel_dy] = run_LG_PKPD(y0,tspan)

%This function numerically simulates the LG PKPD model

global par

%Set up integrator and run
options = odeset('NonNegative',1:length(y0));
[tout,yout] = ode15s(@(t,y) LG_PKPD_RHS(t,y), tspan, y0, options);
yend = yout(end,:);

%Make final call of RHS function to get system state
[dy,kox_eff,~,energy_flux] = LG_PKPD_RHS(1,yend);

rel_dy = dy./yend;


end