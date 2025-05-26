%This script computes the equilibrium concentration of weak acids and
%conjugate base compounds in diets containing lactic and acetic acid.


%% First: how does changing the relative proportion of lactic and acetic
% acid in a mixture with a constant number of carbon atoms. We start with
% 62.5 mM glucose, the 1x CDM concentration. Lactic acid contains three
% carbon atoms, while acetic acid contains two.

clear;clc

%Order is lactic acid, acetic acid
pKa = [3.86;4.76];
H0 = 0;

%Set up for simulation
n = 1000;
f_glu_to_lac = linspace(0,1,n); %Fraction of glucose converted to lactate
C_glu = 0.0625;
C_lac = 2*f_glu_to_lac*C_glu;
C_ace = 3*(1-f_glu_to_lac)*C_glu;

%Set up storage vectors
C_lac_anion = zeros(size(C_lac));
C_ace_anion = zeros(size(C_lac));

%Loop through concentrations
for i = 1:n

    %Initial concentration vectors
    C0 = [C_lac(i);C_ace(i)];
    
    %Solve and store results
    [equi_vec,residuals(i,:),pH(i)] = solve_two_acid_equi(pKa,C0,H0);
    C_lac_anion(i) = equi_vec(1);
    C_ace_anion(i) = equi_vec(2);

end

%Get fraction of compounds as conjugate base
frac_lac_anion = C_lac_anion./C_lac;
frac_ace_anion = C_ace_anion./C_ace;

%Make results table and export
mixed_acid_sim = array2table([f_glu_to_lac',1-frac_ace_anion',1-frac_lac_anion'],...
    'VariableNames',{'frac_glu_to_lac','frac_acetic_acid','frac_lactic_acid'});
writetable(mixed_acid_sim,'mixed_acid_sim.xlsx','WriteVariableNames',true)


%% Look at mixtures just containing lactic acid or acetic acid + glucose

%Order is lactic acid, acetic acid
pKa = [3.86;4.76];
H0 = 0;

%Set up simulations
n = 1000;
f_glu_converted = linspace(0,1,n); %Fraction of glucose converted to acid
C_glu = 0.0625;
C_lac = 2*f_glu_converted*C_glu;
C_ace = 3*f_glu_converted*C_glu;

%Set up storage vectors
C_lac_anion = zeros(size(C_lac));
C_ace_anion = zeros(size(C_lac));

%Loop through concentrations
for i = 1:n

    %Lactic acid simulation
    C01 = [C_lac(i);0];
    [equi_vec,residuals1(i,:),pH1(i)] = solve_two_acid_equi(pKa,C01,H0);
    C_lac_anion(i) = equi_vec(1);

    %Acetic acid simulation
    C02 = [0;C_ace(i)];
    [equi_vec,residuals2(i,:),pH2(i)] = solve_two_acid_equi(pKa,C02,H0);
    C_ace_anion(i) = equi_vec(2);

end

%Get fraction of compounds as conjugate base
frac_lac_anion = C_lac_anion./C_lac;
frac_ace_anion = C_ace_anion./C_ace;

%Generate results table and export
single_acid_sim = array2table([f_glu_to_lac',1-frac_ace_anion',1-frac_lac_anion'],...
    'VariableNames',{'frac_glu_converted','frac_acetic_acid','frac_lactic_acid'});
writetable(mixed_acid_sim,'single_acid_sim.xlsx','WriteVariableNames',true)

