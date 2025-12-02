% VMS Optimization Using Multi-objective Genetic Algorithm
% For 3 Objectives, Energy, Productivity and Water
clear;
clc;

global eval_count;
eval_count = 0;

% ------------------------ Domain Setup ------------------------
N = input("Discretization size (positive integer) = ");
RH = input("Relative Humidity (%) = ");

Pop = 2;  % Populaton size (10*nDV)
Gen = 2;  % Generations (10*log(nDV)

% ------------------------ Optimization Setup -----------------------
% [tads, tvs, tdes, u ads, u vp]
lb = [1000, 1000, 200, 0.1, 0.01];  % Lower bounds: 
ub = [10800, 10800, 400, 1, 10];    % Upper bounds
options = optimoptions('gamultiobj', ...
    'PopulationSize', Pop, ...      % Larger population for more Pareto points
    'MaxGenerations', Gen, ...      % More generations for better convergence
    'ParetoFraction', 0.5, ...      % Retain more Pareto solutions
    'PlotFcn', {@gaplotpareto});%, ...  

% Run Multi-objective Genetic Algorithm
[x_opt, fval] = gamultiobj(@(x) objective_function(x, N, RH), 5, [], [], [], [], lb, ub, [], options);
                                                                                                                                                                  
num_solutions = size(x_opt, 1);
extra_metrics = zeros(num_solutions, 3); 

% Evaluate additional performance metrics for each Pareto solution
for i = 1:num_solutions
    [~, extra_metrics(i, :)] = objective_function(x_opt(i, :), N, RH);
end

% Extract variables for plotting
energy      = fval(:,1) ;                
prodtvy     = -fval(:,2);   
ads_time    = x_opt(:,1);            
vs_time     = x_opt(:,2);   
ds_time     = x_opt(:,3);
v_ads       = x_opt(:,4);         
v_vs        = x_opt(:,5);        
gH2Olost_gCO2    = extra_metrics(:, 1);
gH2Oused_gCO2cap = extra_metrics(:, 2);
CO2_Purity   = extra_metrics(:, 3);

% save("Optimized_P50G60_3Obj_vp10ms_80RHtry_updatedHS.mat",'energy','prodtvy','ads_time','vs_time','ds_time','v_ads','v_vs','gH2Olost_gCO2','gH2Oused_gCO2cap','CO2_Purity','-v7.3' )

% ------------------------ Display Best Solution ------------------------
% Find the best solution based on energy minimization
[~, best_idx] = min(fval(:,1)); % Index of minimum energy

fprintf('\nBest Optimized Solution:\n');
fprintf('Adsorption Time: %.2f sec\n', x_opt(best_idx, 1));
fprintf('Vapor Stripping Time: %.2f sec\n', x_opt(best_idx, 2));
fprintf('Desorption Time: %.2f sec\n', x_opt(best_idx, 3));
fprintf('Velocity Ads: %.2f m/s\n', x_opt(best_idx, 4));
fprintf('Velocity vp : %.2f m/s\n', x_opt(best_idx, 5));
fprintf('Mechanical Energy: %.2f kWh/ton CO2\n', fval(best_idx, 1));
fprintf('CO2 Productivity: %.2f mol/m^3/hr\n', -fval(best_idx, 2));

% ------------------------ Display All Solutions ------------------------
fprintf('\nAll Solutions Found:\n');
fprintf('%10s %10s %10s %5s %5s %15s %20s\n', 't_ad (s)', 't_vp (s)', 't_de (s)', 'u_ads','u_vp','Mech Energy (kWh/ton)', 'CO2 Productivity (mol/m^3/hr)');
for i = 1:size(x_opt, 1)
    fprintf('%10.2f %10.2f %10.2f %5.2f %5.2f %15.2f %20.2f\n', x_opt(i,1), x_opt(i,2), x_opt(i,3),  x_opt(i,4), x_opt(i,5), fval(i,1), -fval(i,2));
end


%% Objective Function
function [objectives, extra_metrics] = objective_function(x, N, RH)
    
    global eval_count;
    eval_count = eval_count + 1;

    t_ad = x(1);  % Adsorption Time
    t_vp = x(2);  % Vapor Stripping Time
    t_de = x(3);  % Desorption Time
    v_ads = x(4);  % Adsorption velocity
    v_vp = x(5);  % Vapor stripping velocity

    % Call the aVMS cycle main script with decision variables
    [MECH_ENERGY, CO2_PRODUCTIVITY, gH2Olost_gCO2cap,gH2Oused_gCO2cap, CO2_Purity] = VMS_MainCycle(N, RH, t_ad, t_vp, t_de, v_ads, v_vp);

    % fprintf('Testing: Adsorption: %.2f, Vapor Stripping: %.2f, Desorption: %.2f, v_ads: %.2f, v_vs: %.2f, Energy: %.2f\n', t_ad, t_vp, t_de, v_ads, v_vp, MECH_ENERGY);
    fprintf('Testing [%d]: Adsorption: %.2f, Vapor Stripping: %.2f, Desorption: %.2f, v_ads: %.2f, v_vs: %.2f, Energy: %.2f\n, Productivity: %.2f\n', ...
            eval_count, t_ad, t_vp, t_de, v_ads, v_vp, MECH_ENERGY, CO2_PRODUCTIVITY);

    % Apply constraint: Mechanical Energy must be ≤ 6000 kWh/ton CO₂
    penalty = 0;
    if MECH_ENERGY > 6000
        penalty = 1e3 * (MECH_ENERGY - 6000)^2;  % Large penalty for exceeding constraint
    end

    % Objectives (to be minimized)
    objectives(1) = MECH_ENERGY + penalty;  % Minimize Mechanical Energy
    objectives(2) = -CO2_PRODUCTIVITY;      % Maximize CO₂ Productivity (negated for minimization)
    objectives(3) = gH2Olost_gCO2cap;

    % Extra metrics
    extra_metrics = [gH2Olost_gCO2cap, gH2Oused_gCO2cap, CO2_Purity];
end


