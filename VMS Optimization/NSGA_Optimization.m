% VMS Optimization Using Multi-objective Genetic Algorithm
% For 3 Objectives: Energy, Productivity and Water
clear;
clc;

% NOTE: In parallel, a global counter will not be consistent across workers.
global eval_count;
eval_count = 0;

% ------------------------ Domain Setup ------------------------
N   = input('Enter the number of spatial nodes, N (e.g., 10): ');
RH  = input('Enter the relative humidity, RH (%) (e.g., 50%): ');

Pop = input('Enter the GA population size (e.g., 2): ');
Gen = input('Enter the number of GA generations (e.g., 3): ');

% ------------------------ Decision Variables -----------------------

% x = [t_cs (s), t_vs (s), t_des (s), v_air (m/s), v_vs (m/s)] 
lb = [1000, 1000, 200, 0.10, 0.01]      ; % Lower bounds
ub = [10800, 10800, 500, 1.00, 5.00]    ; % Upper bounds
x0 = [5000, 10000, 200, 0.50, 1.00]     ; % Initial Guess

% ------------------------ Parallel Setup ------------------------
pool = gcp('nocreate');
if isempty(pool)
    parpool;
end

options = optimoptions('gamultiobj', ...
    'PopulationSize', Pop, ...
    'MaxGenerations', Gen, ...
    'ParetoFraction', 0.5, ...
    'InitialPopulationMatrix', x0, ...
    'UseParallel', true, ...
    'Display', 'iter', ...
    'PlotFcn', []);

% Run Multi-objective Genetic Algorithm
[x_opt, fval] = gamultiobj(@(x) objective_function(x, N, RH), ...
    5, [], [], [], [], lb, ub, options);

num_solutions = size(x_opt, 1);
extra_metrics = zeros(num_solutions, 3);

% Evaluate additional performance metrics for each Pareto solution
parfor i = 1:num_solutions
    [~, extra_metrics(i, :)] = objective_function(x_opt(i, :), N, RH);
end

% Extract variables for plotting
energy     = fval(:,1)      ; % (kWh/tCO2)
prodtvy    = -fval(:,2)     ; % (mol/m3/hr)
cs_time    = x_opt(:,1)     ; % (s)
vs_time    = x_opt(:,2)     ; % (s)
des_time   = x_opt(:,3)     ; % (s)
vair       = x_opt(:,4)     ; % (m/s)
vvs        = x_opt(:,5)     ; % (m/s)
gH2Olost_gCO2    = extra_metrics(:, 1); % g/g
gH2Oused_gCO2cap = extra_metrics(:, 2); % g/g
CO2_Purity       = extra_metrics(:, 3); % g/g

save("Optimized_P50G30_3Obj_50RH_Mk_EnCorrection_fullyloaded.mat", ...
    'energy','prodtvy','cs_time','vs_time','des_time','vair','vvs', ...
    'gH2Olost_gCO2','gH2Oused_gCO2cap','CO2_Purity','-v7.3')


% ------------------------ Display Best Solution ------------------------
[~, best_idx] = min(fval(:,1)); % Index of minimum energy

fprintf('\nBest Optimized Solution  :\n');
fprintf('Sorption Time              : %.2f sec\n', x_opt(best_idx, 1));
fprintf('Vapor Stripping Time       : %.2f sec\n', x_opt(best_idx, 2));
fprintf('Desorption Time            : %.2f sec\n', x_opt(best_idx, 3));
fprintf('Air velocity               : %.2f m/s\n', x_opt(best_idx, 4));
fprintf('VS velocity                : %.2f m/s\n', x_opt(best_idx, 5));
fprintf('Mechanical Energy          : %.2f kWh/ton CO2\n', fval(best_idx, 1));
fprintf('CO2 Productivity           : %.2f mol/m^3/hr\n', -fval(best_idx, 2));

% ------------------------ Display All Solutions ------------------------
fprintf('\nAll Solutions Found:\n');
fprintf('%10s %10s %10s %5s %5s %15s %20s %20s\n', ...
    't_cs (s)', 't_vs (s)', 't_des (s)', 'v_air','v_vs', ...
    'Mech Energy (kWh/ton)', 'CO2 Productivity (mol/m^3/hr)', 'gH2Olost/gCO2');

for i = 1:size(x_opt, 1)
    fprintf('%10.2f %10.2f %10.2f %5.2f %5.2f %15.2f %20.2f %20.2f\n', ...
        x_opt(i,1), x_opt(i,2), x_opt(i,3), x_opt(i,4), x_opt(i,5), ...
        fval(i,1), -fval(i,2), fval(i,3));
end

%% Objective Function (with CSS penalty + max cycles passed to VMS_MainCycle)
function [objectives, extra_metrics] = objective_function(x, N, RH)
    global eval_count;
    eval_count = eval_count + 1;

    t_cs  = x(1);
    t_vs  = x(2);
    t_des = x(3);
    vair  = x(4);
    vvs   = x(5);

    % ---------- CSS settings ---------
    BIG = 1e9;                       % <-- domination penalty for non-CSS cases

    % Default outputs (in case of failure)
    objectives = [BIG, BIG, BIG];
    extra_metrics = [NaN, NaN, NaN];

    try        
        [MECH_ENERGY, CO2_PRODUCTIVITY, gH2Olost_gCO2cap, ...
         gH2Oused_gCO2cap, CO2_Purity, CSS_ok,i_end] = ...
            VMS_MainCycle(N, RH, t_cs, t_vs, t_des, vair, vvs);

        % ---------- DISCARD non-CSS cases ----------
        % If CSS is not achieved within maxCycles (i_end==maxCycles), apply huge penalty
        if ~CSS_ok
            objectives = [BIG, BIG, BIG];
            extra_metrics = [NaN, NaN, NaN];
            return;
        end

        % ---------- Energy constraint penalty ----------
        penalty = 0;
        if MECH_ENERGY > 20000
            penalty = 1e3 * (MECH_ENERGY - 20000)^2;
        end

        % ---------- Objectives  ----------
        objectives(1) = MECH_ENERGY + penalty;
        objectives(2) = -CO2_PRODUCTIVITY;
        objectives(3) = gH2Olost_gCO2cap;

        % ---------- Extra metrics ----------
        extra_metrics = [gH2Olost_gCO2cap, gH2Oused_gCO2cap, CO2_Purity];

    catch ME
        % Any crash/numerical issue becomes a discarded (dominated) point
        fprintf('Eval [%d] failed: %s\n', eval_count, ME.message);
        objectives = [BIG, BIG, BIG];
        extra_metrics = [NaN, NaN, NaN];
    end
end
