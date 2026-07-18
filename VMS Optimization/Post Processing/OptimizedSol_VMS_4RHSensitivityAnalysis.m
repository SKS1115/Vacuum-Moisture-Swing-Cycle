%% Clear workspace and figures
clear; clc; close all;

%% Load Excel data (4 RH cases)
files = { ...
    'Optimized_P50G30_3Obj_20RH_Mk_EnCorrection_fullyloaded.xlsx', ...
    'Optimized_3obj_P50G30_50RH_MK_EnCorrection_fullyloaded.xlsx', ...
    'Optimized_P50G30_3Obj_65RH_Mk_EnCorrection_fullyloaded.xlsx', ...
    'Optimized_3obj_P50G30_80RH_MK_EnCorrection.xlsx' ...
    };

RH_labels = {'20% RH', '50% RH', '65% RH', '80% RH'};

% 4 distinct markers
symbols = {'o', 's', 'd', '^'};

% Read data
N = numel(files);

data = cell(1,N);
for i = 1:N
    data{i} = readtable(files{i});
end

%% Extract variables
Energy = cell(1,N);
Productivity = cell(1,N);
H2O_Lost_CO2 = cell(1,N);

for i = 1:N
    Energy{i} = data{i}.Energy;
    Productivity{i} = data{i}.Productivity;
    H2O_Lost_CO2{i} = data{i}.("gH2OLost_gCO2"); % adjust if header differs
end

%% === Figure: Fixed 0–4.5 color scale for water lost ===
figure(1); clf;

axMain = axes('Position',[0.10 0.12 0.75 0.80]); % room on right for colorbars
hold(axMain,'on');

% % Define color gradients (light → dark)
% brickRed  = [linspace(1.0,0.6,256)', linspace(0.8,0.1,256)', linspace(0.8,0.1,256)'];   % reddish
% purpleMap = [linspace(0.9,0.4,256)', linspace(0.75,0.1,256)', linspace(1.0,0.6,256)']; % purple
% blueMap   = [linspace(0.7,0.0,256)', linspace(0.9,0.2,256)', linspace(1.0,0.3,256)'];   % blue
% greenMap  = [linspace(0.6,0.1,256)', linspace(0.9,0.4,256)', linspace(0.6,0.1,256)'];   % green/teal
% 
% % Use 4 maps only
% cmap_list = {brickRed, purpleMap, blueMap, greenMap};

nC = 256;
t  = linspace(0,1,nC)';   % 0 = low water loss, 1 = high water loss

% Make the low end whiter by mixing with white more strongly at small t
% (increase "whiten_strength" to get even whiter lows)
whiten_strength = 0.8;                 % 0–1 (0 = no whitening, 1 = very white)
w = whiten_strength*(1 - t).^1.2;       % whitening weight (more at low end)

% Base (your original) gradients
brickBase  = [linspace(1.0,0.6,nC)', linspace(0.8,0.1,nC)', linspace(0.8,0.1,nC)'];
purpleBase = [linspace(0.9,0.4,nC)', linspace(0.75,0.1,nC)', linspace(1.0,0.6,nC)'];
blueBase   = [linspace(0.7,0.0,nC)', linspace(0.9,0.2,nC)', linspace(1.0,0.3,nC)'];
greenBase  = [linspace(0.6,0.1,nC)', linspace(0.9,0.4,nC)', linspace(0.6,0.1,nC)'];

% Blend with white near the low end
white = ones(nC,3);
brickRed  = (1-w).*brickBase  + w.*white;
purpleMap = (1-w).*purpleBase + w.*white;
blueMap   = (1-w).*blueBase   + w.*white;
greenMap  = (1-w).*greenBase  + w.*white;

% Use 4 maps
cmap_list = {brickRed, purpleMap, blueMap, greenMap};

% Fixed color scale for water lost
water_min = 0;
water_max = 4;

scatterHandles = gobjects(1,N);

for i = 1:N
    cvals = H2O_Lost_CO2{i};
    cmap  = cmap_list{i};

    % Map data to colors using fixed scale
    color_idx = round(rescale(cvals, 1, size(cmap,1), ...
        'InputMin', water_min, 'InputMax', water_max));
    color_idx = max(min(color_idx, size(cmap,1)), 1); % clamp
    colors = cmap(color_idx, :);

    % scatterHandles(i) = scatter3(axMain, Productivity{i}, Energy{i}, H2O_Lost_CO2{i}, ...
    %     150, colors, symbols{i}, 'filled', ...
    %     'DisplayName', RH_labels{i});

    scatterHandles(i) = scatter3(axMain, Productivity{i}, Energy{i}, H2O_Lost_CO2{i}, ...
    150, colors, symbols{i}, 'filled', ...
    'MarkerEdgeColor', 'k', ...     % black outline
    'LineWidth', 0.5, ...            % thickness of outline
    'DisplayName', RH_labels{i});
end

xlabel(axMain,'Productivity [kg_{CO_2}/kg_s/day]','FontWeight','bold','FontSize',20);
ylabel(axMain,'Electrical Energy [MJ/kg CO_2]','FontWeight','bold','FontSize',20);
zlabel(axMain,'gH_2O_{lost}/gCO_2','FontWeight','bold','FontSize',20);

set(axMain,'FontSize',15,'FontWeight','bold','LineWidth',2);
grid(axMain,'on'); box(axMain,'on');
xlim(axMain,[0 1]); ylim(axMain,[0 30]);

legend(axMain,'Location','best');

hold(axMain,'off');

%% Add mini colorbars for each RH case

for i = 1:4
    % Create small invisible axes for each colorbar
    ax_cb = axes('Position', [0.81, 0.75-(i-1)*0.18, 0.02, 0.12]);
    cmap = cmap_list{i};
    imagesc(linspace(water_min, water_max, size(cmap,1))');
    colormap(ax_cb, cmap);
    ax_cb.YDir = 'normal';
    ax_cb.XTick = [];
    ax_cb.Box = 'on';
    ax_cb.YTick = linspace(1, size(cmap,1), 3);
    yticklabels(ax_cb, sprintfc('%.1f', linspace(water_min, water_max, 3)));
    ylabel(ax_cb, RH_labels{i}, 'FontWeight', 'bold', 'FontSize', 12);
end

