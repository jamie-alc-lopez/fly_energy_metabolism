%% This code analyzes the growth rate predictions from iDrosophila when
% supplied fixed fluxes of lactic acid or glucose. The code for iDrosophila
% configuration is based on codes graciously supplied by Kiran Patil and
% Müberra Fatma.

% Note that this code requires the COBRA toolbox and installation of the
% Gurobi solver.

clear;clc

%Initialize the COBRA toolbox
initCobraToolbox

%Load the iDrosophila1 model struct
model = load('iDrosophila_materials/iDrosophila1.mat');
model = model.iDrosophila1;

% Setting the boundaries of ATPM reaction (non-growth-associated
% maintenance from Schönborn et al., 2019)
model.lb(strcmp('HMR_3964', model.rxns)) = 8.55;
model.ub(strcmp('HMR_3964', model.rxns)) = 8.55;

% Load in Holidic diet from Schönborn et al., 2019
[~, HD_medium, ~] = xlsread('iDrosophila_materials/Holidic_diet.xlsx', ...
    'HD_growth_medium','D5:D51');

%Find exchange reactions corresponding to Holidic diet inputs
HD_ind = arrayfun(@(counter) find(strcmp(HD_medium{counter}, model.equations)),...
    1:length(HD_medium))';

%Find exchange indices
[sel_exc, sel_upt] = findExcRxns(model);
exc_ind = find(sel_exc == 1);
upt_ind = find(sel_upt == 1);
exc_rxns = model.equations(sel_exc);
sec_ind = setdiff(exc_ind, upt_ind);

%Blocking all uptake reactions (but free secretion)
model.lb(exc_ind) = -1000;
model.ub(exc_ind) = 0;
model.lb(sec_ind) = 0;
model.ub(sec_ind) = 1000;
model.lb(strcmp(' <=> O2 ', model.equations)) = 0;
model.ub(strcmp(' <=> O2 ', model.equations)) = 1000;

%Initial Holidic diet max sucrose uptake
Vsuc = 2.21240473;

%Most other metabolites (e.g., aminoacids): Vsuc*1/10
model.ub(HD_ind) = Vsuc/10;

%Identify vitamin uptake and set to Vsuc*1/100
vit_list = {' <=> pantothenate ', ' <=> nicotinate ', ' <=> riboflavin ', ' <=> folate ',...
    ' <=> thiamin ', ' <=> biotin ', ' <=> pyridoxine ', ' <=> gamma-tocopherol ',...
    ' <=> retinoate ', ' <=> alpha-tocopherol ', ' <=> aquacob(III)alamin '}';
vit_ind = arrayfun(@(counter) find(strcmp(vit_list(counter), model.equations)),...
    1:length(vit_list));
model.ub(vit_ind) = Vsuc*1/100;

%Set maximum oxygen uptake
model.ub(strcmp(' <=> O2 ', model.equations)) = 24;

%Allow unlimited uptake of H2O, O2, and salt ions
free_uptake_list = {' <=> Ca2+ '; ' <=> Fe2+ '; ' <=> K+ '; ' <=> Na+ '; ' <=> H2O ';...
    ' <=> zinc '; ' <=> Pi '; ' <=> sulfate '};
free_uptake_ind = arrayfun(@(counter) find(strcmp(free_uptake_list(counter),model.equations)),...
    1:length(free_uptake_list));
model.ub(free_uptake_ind) = 1000;

%Departure from Holidic diet: zero out the sucrose exchange reaction, to be
%replaced by the glc and lac reactions
suc_transp_ind = strcmp(' <=> sucrose ', model.equations);
model.lb(suc_transp_ind) = 0;
model.ub(suc_transp_ind) = 0;

%Identify the lactate and glucose exchange reactions
lac_transp_ind = strcmp(' <=> L-lactate ', model.equations);
glc_transp_ind = strcmp(' <=> glucose ', model.equations);


%% Make model struct to feed to Gurobi, test with default HD configuration

%Make new model struct
modelD.S = sparse(model.S);
modelD.A = modelD.S;
modelD.ub = model.ub;
modelD.lb = model.lb;
modelD.rhs = zeros(length(model.mets),1);
modelD.mets = model.mets;
modelD.metNames = model.metNames;
modelD.rxns = model.rxns;
modelD.equations = model.equations;
modelD.sense = '=';
modelD.vtype = 'C';
modelD.genes = model.genes;
modelD.obj = zeros(length(model.rxns),1);
biom = find(strcmp('Biomass_formation', model.rxns));
modelD.obj(biom) = -1;
modelD.rules = model.rules;
modelD.c = model.c;

%Set max sucrose uptake to Vsuc
modelD.ub(suc_transp_ind) = Vsuc;

%Solve, growth under these conditions should be 0.7738
HD_V = gurobi(modelD);
HD_growth = HD_V.x(biom);

%% Now loop through and estimate growth rate for different input levels of
%glucose and lactate

%Nutrient levels
n1 = 2;
nutr_vec = [1 2];

%Fraction of glucose
n2 = 50;
frac_vec = linspace(0,1,n2);

%Zero out all major carbon inputs
modelD.ub(suc_transp_ind) = 0;
modelD.lb(lac_transp_ind) = 0;
modelD.lb(glc_transp_ind) = 0;

%Loop through nutrient levels
for i = 1:n1

    %Get nutrient level i
    nutr_level = nutr_vec(i);

    %Loop through diet compositions
    for j = 1:n2

        %Get the jth glucose fraction
        frac_j = frac_vec(j);

        %Set maximum input flux. Input carbon flux is set based on the sucrose
        %flux in the Holidic diet. Sucrose is a dissacharaide, so glucose gets
        %double the flux and lactate gets quadruple the flux
        modelD.ub(lac_transp_ind) = nutr_level*(1-frac_j)*4*Vsuc;
        modelD.ub(glc_transp_ind) = nutr_level*frac_j*2*Vsuc;

        %Optional: fix lower bound flux as well
        %modelD.lb(lac_transp_ind) =  nutr_level*(1-frac_j)*4*Vsuc;
        %modelD.lb(glc_transp_ind) = nutr_level*frac_j*2*Vsuc;

        %Solve model
        V{i,j} = gurobi(modelD);

        %Save growth rate
        Vgrowth(i,j) = V{i,j}.x(biom);

    end

end

% %Optional: make figure to visualize behavior
% figure
% hold on
% plot(frac_vec,Vgrowth(1,:),'k-','LineWidth',1.5)
% plot(frac_vec,Vgrowth(2,:),'r-','LineWidth',1.5)
% ylim([0,1.3])
% xlabel('Fraction of glucose (balance lactic acid)')
% ylabel('iDrosophila predicted growth (1/h)')
% set(gca,'FontSize',15)
% legend({'1x','2x'},'FontSize',15)
% name_vector = ['iDrosophila_glc_lac_plot.pdf'];
% exportgraphics(gcf,name_vector,'ContentType','vector')


%Make results table and export
out_table = table(frac_vec',Vgrowth(1,:)',Vgrowth(2,:)','VariableNames',...
    {'frac glc','1x growth','2x growth'});
writetable(out_table,'iDrosophila_glc_lac_results.xlsx');
