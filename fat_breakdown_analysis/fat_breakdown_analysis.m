%% This script analyze the fat breakdown assay data. It does this via three
% fitting analyses. The first fit uses a piecewise linear model. It
% identifies interval in which half utilization is crossed and fitting a
% linear model to that interval. It then estimates the 50% crossing time 
% and propagates the uncertainty. It then performs the same t50 estimation 
% using a linear fit of the full timecourse. It then fits an exponential 
% model
%This script also performs a t-test between fat levels at days 2 and 5,
%comparing all conditions to the glucose value

% Note: this script assumes the 50% point is crossed only once in data.

clear;clc

% Import fat breakdown data, identify conditions for fits
data = readtable('fat_breakdown_data.xlsx');
conditions = {'stv','glc','lac','glc_lac'};
fed_ind = strcmp(data.condition,'fed');
test_days = [2 5];
ref_data = data(strcmp(data.condition,'glc'),:);

yc = 50; %percentage of fat remaining to focus on 
x0 = 1.5; %initial point for root finding

%Loop through different conditions
for i = 1:length(conditions)

    %Get data corresponding to condition
    condition_i = conditions{i};
    data_i = data(strcmp(data.condition,condition_i)|fed_ind,:);

    %Compute t-tests w.r.t. glucose
    for j = 1:length(test_days)
        ref_data_j = ref_data(ref_data.time==test_days(j),:);
        test_data_ij = data_i(data_i.time==test_days(j),:);

        if ~isempty(test_data_ij)
            [~,p(i,j)] = ttest2(ref_data_j.fat_val,test_data_ij.fat_val,...
                'Vartype','unequal','Tail','both');
        else
            p(i,j) = NaN;
        end

    end

    
    %Perform the piece linear fit analysis ---------------------

    %Average data
    G = findgroups(data_i.condition,data_i.time);
    t = splitapply(@mean,data_i.time,G);
    y = splitapply(@mean,data_i.fat_val,G);

    %Fit a linear interpolation model to the average
    mean_interp_fit = fit(t,y,'linearinterp');
    [init_tc,fval] = fsolve(@(x) (mean_interp_fit(x)-yc).^2,x0);

    %Use interpolation to find the interval containing the 50% point
    above_ind = find(t>=init_tc);
    above_ind = above_ind(1);
    time_1 = t(above_ind);
    below_ind = find(t < init_tc);
    below_ind = below_ind(end);
    time_2 = t(below_ind);

    %Get data from the interval containing the 50% point
    interval_i = data_i(data_i.time == time_1 | data_i.time == time_2,:);

    %Fit and solve interval linear model
    interval_fit = fit(interval_i.time,interval_i.fat_val,'poly1');
    p1_interval = interval_fit.p1;
    p2_interval = interval_fit.p2;
    tc_interval_linear(i) = (yc - p2_interval)/p1_interval;

    %Extract errors
    c_ints = confint(interval_fit);
    sigma_p1_interval = (c_ints(2,1) - c_ints(1,1))/2;
    sigma_p2_interval = (c_ints(2,2) - c_ints(1,2))/2;

    %Compute partial derivatives and propagate error
    diff_tc_p1_interval = -1/p1_interval;
    diff_tc_p2_interval = (p2_interval - yc)./(p1_interval^2);
    sigma_tc_interval_linear(i) = sqrt(sigma_p1_interval^2*diff_tc_p1_interval^2 ...
        + sigma_p2_interval^2*diff_tc_p2_interval^2);


    %Now, do a full linear fit -------------------------
    lower = [-inf,100]; %Fix intercept at 100
    upper = [inf,100];
    [full_linear_fit,full_linear_gof] = ...
        fit(data_i.time,data_i.fat_val,'poly1','Lower',lower,'Upper',upper);
    full_linear_adj_rsquared(i) = full_linear_gof.adjrsquare;

    full_linear_rate(i) = full_linear_fit.p1;
    p2_full = full_linear_fit.p2;    %Average data
    G = findgroups(data_i.condition,data_i.time);
    t = splitapply(@mean,data_i.time,G);
    y = splitapply(@mean,data_i.fat_val,G);
    tc_full_linear(i) = (yc - p2_full)/full_linear_rate(i);

    %Extract errors
    c_ints = confint(full_linear_fit);
    sigma_p1_full(i) = (c_ints(2,1) - c_ints(1,1))/2;
    sigma_p2_full = 0;

    %Compute partial derivatives and propagate error
    diff_tc_p1_full = -1/full_linear_rate(i);
    diff_tc_p2_full = (p2_full - yc)./(full_linear_rate(i)^2);
    sigma_tc_full_linear(i) = ...
        sqrt(sigma_p1_full(i)^2*diff_tc_p1_full^2 + sigma_p2_full^2*diff_tc_p2_full^2);


    %Perform exponential fit -----------------------
    lower = [100,-inf];
    upper = [100,inf];
    startpoint = [100,-0.4];
    [exp_fit,exp_gof] = ...
        fit(data_i.time,data_i.fat_val,'exp1','Lower',lower,'Upper',upper,...
        'StartPoint',startpoint);
    
    exp_rate(i) = exp_fit.b;
    exp_adj_rsquared(i) = exp_gof.adjrsquare;

end

%Make table and export
fat_breakdown_stats = ...
    table(conditions',p(:,1),p(:,2),full_linear_rate',sigma_p1_full',...
    full_linear_adj_rsquared',tc_full_linear',sigma_tc_full_linear',...
    tc_interval_linear',sigma_tc_interval_linear',exp_rate',...
    exp_adj_rsquared',...
    'VariableNames',{'condition','p_day2','p_day5','linear_rate','sigma_linear_rate',...
    'full_linear_adj_r2','t50_full_linear','sigma_t50_full',...
    't50_interval','sigma_t50_interval','exp_rate','exp_adj_r2'});

writetable(fat_breakdown_stats,['fat_breakdown_stats.xlsx'])