% 

%% Clear workspace and figures
clear; clc; close all;

%% Load Excel data
files = {
    'Optimized_3obj_P50G30_20RH_Sk_EnCorrection.xlsx', ...
    'Optimized_3obj_P50G30_20RH_Mk_EnCorrection.xlsx', ...
    'Optimized_3obj_P50G30_20RH_Fk_EnCorrection.xlsx'
};
kinetic_labels = {'slow k_i', 'measured k_i', 'fast k_i'};
colorbar_energy = {'Energy','Energy','Energy'};
colorbar_water = {'gH_2O_{lost}/gCO_2','gH_2O_{lost}/gCO_2','gH_2O_{lost}/gCO_2'};
symbols = {'o', 's', 'd'};  % Circle, square, diamond for each case

data = cell(1,3);
for i = 1:3
    data{i} = readtable(files{i});
end

%% Extract variables
Energy = cell(1,3);
Productivity = cell(1,3);
H2O_Lost_CO2 = cell(1,3);

for i = 1:3
    Energy{i} = data{i}.Energy;
    Productivity{i} = data{i}.Productivity;
    H2O_Lost_CO2{i} = data{i}.("gH2OLost_gCO2"); % Adjust if header differs
end

%% === Figure 1: Orange, Dark Green, and Dark Blue color gradients ===
figure(1); clf; hold on;

% Define custom colormaps (light → dark)
orangeMap = [linspace(1.0, 0.8, 256)', linspace(0.8, 0.3, 256)', linspace(0.5, 0.0, 256)'];  % Light orange → deep orange
darkGreenMap = [linspace(0.7, 0.0, 256)', linspace(1.0, 0.4, 256)', linspace(0.7, 0.0, 256)']; % Light green → dark green
darkBlueMap = [linspace(0.7, 0.0, 256)', linspace(0.85, 0.1, 256)', linspace(1.0, 0.3, 256)']; % Light blue → dark blue

cmap_list = {orangeMap, darkGreenMap, darkBlueMap};

% Global scaling for water lost: 0–7 (light → dark)
water_min = 0;
water_max = max(cellfun(@max, H2O_Lost_CO2));

scatterHandles = gobjects(1,3);
for i = 1:3
    % Normalize colors using fixed 0–7 range
    cvals = H2O_Lost_CO2{i};
    cmap = cmap_list{i};
    color_idx = round(rescale(cvals, 1, size(cmap,1), 'InputMin', water_min, 'InputMax', water_max));
    color_idx = max(min(color_idx, size(cmap,1)), 1); % Clamp indices
    colors = cmap(color_idx, :);

    scatterHandles(i) = scatter3(Productivity{i}, Energy{i}, H2O_Lost_CO2{i}, ...
        150, colors, symbols{i}, 'filled', 'DisplayName', kinetic_labels{i});
end

xlabel('Productivity [kg {CO_2}/kg_s/day]', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('Electrical Energy [MJ/kg CO_2]', 'FontWeight', 'bold', 'FontSize', 20);
zlabel('g{H_2O}_{lost}/g{CO_2}', 'FontWeight', 'bold', 'FontSize', 20);
set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'LineWidth', 2);
grid on; box on; xlim([0 1.2]); ylim([0 5]);
legend('Location', 'best');
% title('Figure 1: Color Gradients (0–7 g{H_2O}_{lost}/g{CO_2})', 'FontWeight', 'bold', 'FontSize', 18);

%% Add colorbars for each case (same 0–7 scale)
for i = 1:3
    ax_cb = axes('Position', [0.88, 0.75-(i-1)*0.22, 0.02, 0.15]);
    cmap = cmap_list{i};
    imagesc(linspace(water_min, water_max, size(cmap,1))');
    colormap(ax_cb, cmap);
    ax_cb.YDir = 'normal';
    ax_cb.XTick = [];
    ax_cb.Box = 'on';
    ax_cb.YTick = linspace(1, size(cmap,1), 3);
    yticklabels(ax_cb, sprintfc('%.1f', linspace(water_min, water_max, 3)));
    ylabel(ax_cb, colorbar_water{i}, 'FontWeight', 'bold', 'FontSize', 12);
end

hold off;

%% Figure 2: Water loss vs Productivity for all kinetic cases

figure(2); clf; hold on;

% Define custom colormaps (light → dark)
orangeMap = [linspace(1.0, 0.8, 256)', linspace(0.8, 0.3, 256)', linspace(0.5, 0.0, 256)'];  % Light orange → deep orange
darkGreenMap = [linspace(0.7, 0.0, 256)', linspace(1.0, 0.4, 256)', linspace(0.7, 0.0, 256)']; % Light green → dark green
darkBlueMap = [linspace(0.7, 0.0, 256)', linspace(0.85, 0.1, 256)', linspace(1.0, 0.3, 256)']; % Light blue → dark blue

cmap_list = {orangeMap, darkGreenMap, darkBlueMap};

% Global scaling for energy: adjust based on your data range
energy_min = 0.8;%min(cellfun(@min, Energy));
energy_max = 30;%max(cellfun(@max, Energy));

scatterHandles = gobjects(1,3);
for i = 1:3
    % Normalize colors using fixed energy range
    cvals = Energy{i};
    cmap = cmap_list{i};
    color_idx = round(rescale(cvals, 1, size(cmap,1), 'InputMin', energy_min, 'InputMax', energy_max));
    color_idx = max(min(color_idx, size(cmap,1)), 1); % Clamp indices
    colors = cmap(color_idx, :);

    scatterHandles(i) = scatter(Productivity{i}, H2O_Lost_CO2{i}, ...
        150, colors, symbols{i}, 'filled', 'DisplayName', kinetic_labels{i});
end

xlabel('Productivity [kg {CO_2}/kg_s/day]', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('g{H_2O}_{lost}/g{CO_2}', 'FontWeight', 'bold', 'FontSize', 20);
xlim([0 1.2]); ylim([0.5 5]);
set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'LineWidth', 2);
grid on; box on;
legend('Location', 'best');

%% Add colorbars for each case (based on energy scale)
for i = 1:3
    ax_cb = axes('Position', [0.88, 0.75-(i-1)*0.22, 0.02, 0.15]);
    cmap = cmap_list{i};
    imagesc(linspace(energy_min, energy_max, size(cmap,1))');
    colormap(ax_cb, cmap);
    ax_cb.YDir = 'normal';
    ax_cb.XTick = [];
    ax_cb.Box = 'on';
    ax_cb.YTick = linspace(1, size(cmap,1), 3);
    yticklabels(ax_cb, sprintfc('%.1f', linspace(energy_min, energy_max, 3)));
    ylabel(ax_cb, colorbar_energy{i}, 'FontWeight', 'bold', 'FontSize', 12);
end

hold off;


% %% === Figure 2: Energy vs Adsorption Time for all kinetic cases ===
% figure(2); clf; hold on;
% 
% for i = 1:3
%     % Color represents Adsorption Velocity
%     cvals = data{i}.V_Ads;  
%     cmap = cmap_list{i};
%     color_idx = round(rescale(cvals, 1, size(cmap,1)));
%     color_idx = max(min(color_idx, size(cmap,1)), 1);
%     colors = cmap(color_idx, :);
% 
%     scatter(data{i}.AdsTime, Energy{i}, 150, colors, symbols{i}, 'filled', 'DisplayName', kinetic_labels{i});
% end
% 
% xlabel('Adsorption Time [s]', 'FontWeight', 'bold', 'FontSize', 20);
% ylabel('Electrical Energy [MJ/kg CO_2]', 'FontWeight', 'bold', 'FontSize', 20);
% set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'LineWidth', 2);
% grid on; box on;
% legend('Location', 'best');
% title('Figure 2: Energy vs Adsorption Time', 'FontWeight', 'bold', 'FontSize', 18);
% 
% % === Individual colorbars for each kinetic case (0–1 range) ===
% for i = 1:3
%     colormap(gca, cmap_list{i});
%     cb = colorbar('Position', [0.93 + (i-1)*0.03, 0.15, 0.02, 0.7]); % shift each bar right
%     caxis([0 1]);
%     cb.Label.String = sprintf('V_{Ads} [m/s] (%s)', kinetic_labels{i});
%     cb.Label.FontSize = 13;
%     cb.FontWeight = 'bold';
%     cb.FontSize = 11;
% end
% 
% hold off;
% 
% 
% %% === Figure 3: Productivity vs Vapor Stripping Time for all kinetic cases ===
% figure(3); clf; hold on;
% 
% for i = 1:3
%     % Color represents Vapor Stripping Velocity
%     cvals = data{i}.V_VS;  
%     cmap = cmap_list{i};
%     color_idx = round(rescale(cvals, 1, size(cmap,1)));
%     color_idx = max(min(color_idx, size(cmap,1)), 1);
%     colors = cmap(color_idx, :);
% 
%     scatter(data{i}.VSTime, Productivity{i}, 150, colors, symbols{i}, 'filled', 'DisplayName', kinetic_labels{i});
% end
% 
% xlabel('Vapor Stripping Time [s]', 'FontWeight', 'bold', 'FontSize', 20);
% ylabel('Productivity [kg_{CO_2}/kg_s/day]', 'FontWeight', 'bold', 'FontSize', 20);
% set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'LineWidth', 2);
% grid on; box on;
% legend('Location', 'best');
% title('Figure 3: Productivity vs Vapor Stripping Time', 'FontWeight', 'bold', 'FontSize', 18);
% 
% % === Individual colorbars for each kinetic case (0–10 range) ===
% for i = 1:3
%     colormap(gca, cmap_list{i});
%     cb = colorbar('Position', [0.93 + (i-1)*0.03, 0.15, 0.02, 0.7]); % spaced colorbars
%     caxis([0 10]);
%     cb.Label.String = sprintf('V_{VS} [m/s] (%s)', kinetic_labels{i});
%     cb.Label.FontSize = 13;
%     cb.FontWeight = 'bold';
%     cb.FontSize = 11;
% end
% 
% hold off;
