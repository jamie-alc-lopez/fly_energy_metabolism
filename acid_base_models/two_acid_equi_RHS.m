function RHS = two_acid_equi_RHS(x,Ka,C0,H0)

%This function computes the RHS of the acid-base equilibrium expression for
%a mixture of two weak acids

RHS = zeros(size(Ka));

RHS(1) = x(1)*(H0 + x(1) + x(2)) - Ka(1)*(C0(1) - x(1));
RHS(2) = x(2)*(H0 + x(1) + x(2)) - Ka(2)*(C0(2) - x(2));

end