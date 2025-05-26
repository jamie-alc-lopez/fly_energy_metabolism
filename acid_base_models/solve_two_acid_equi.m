function [equi_vec,residuals,pH] = solve_two_acid_equi(pKa,C0,H0)

%This function solves for the acid/base equilibrium of a mixture of two
%acids

%Get Ka from pKa
Ka = 10.^(-pKa);

%Define RHS function
RHS_fun = @(x) two_acid_equi_RHS(x,Ka,C0,H0)./Ka;

%Set initial conditions and solver options
x0 = 0.5*C0;
options = optimoptions('fsolve','FunctionTolerance',1e-10);

%Solve 
equi_vec = fsolve(RHS_fun,x0,options);
residuals = RHS_fun(equi_vec).*Ka;

%Compute pH
pH = -log10(H0 + sum(equi_vec));

end