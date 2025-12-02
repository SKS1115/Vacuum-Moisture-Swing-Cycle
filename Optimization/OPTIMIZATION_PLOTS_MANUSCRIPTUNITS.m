


% load Optimized_P50G10_NoEnergyCorr_3Obj_vp10ms_20RH.mat


% Define blue color
blueColor = [0 0.4470 0.7410];

%% Figure 1: Energy vs Productivity
% Convert units for manuscript
energy  = energy * (3.6/1000)           ; % from kWh/tonCO2 to MJ/kgCO2
prodtvy = prodtvy *((24*0.044)/(1060))  ; % from mol/m3/hr to kgCO2/kgs/day
ads_time= ads_time/60                   ; % convert s to min
vs_time = vs_time/60                    ;   % convert s to min
ds_time = ds_time/60                    ; % convert s to min

% Combine all variables into a table
% T = table(energy(:), prodtvy(:), ads_time(:), vs_time(:), ds_time(:), ...
%           v_ads(:), v_vs(:), gH2Olost_gCO2(:), gH2Oused_gCO2cap(:), CO2_Purity(:), ...
%           'VariableNames', {'Energy', 'Productivity', 'AdsTime', 'VSTime', 'DSTime', ...
%                             'V_Ads', 'V_VS', 'gH2OLost_gCO2', 'gH2OUsed_gCO2Cap', 'CO2 Purity'});
% writetable(T, 'Optimized_3obj_P50G80_adjustedHS_8-29-25.xlsx');

figure(1); clf;
scatter3(prodtvy, energy, gH2Olost_gCO2, 150, gH2Olost_gCO2, 'filled'); colormap(sky); colorbar;
xlabel('Productivity [kg {CO_2}/kg_s/day]', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('Electrical Energy [MJ/kg CO_2]', 'FontWeight', 'bold', 'FontSize', 20);
zlabel('g{H_2O}_{lost}/g{CO_2}', 'FontWeight', 'bold', 'FontSize', 20);
c = colorbar; c.Label.String = 'g{H_2O}_{lost}/g{CO_2}'; c.Label.FontSize = 15; c.Label.FontWeight = 'bold'; 
set(gca, 'FontSize', 15, 'FontWeight', 'bold'); % Make axis markers bold
grid on; set(gca, 'LineWidth', 2); hold off;
% ylim([0 6]); xlim([0.05 0.25]); zlim([2 7])

figure(2); clf;
scatter3(prodtvy, energy, gH2Olost_gCO2, 150, energy, 'filled'); colormap(flipud(autumn)); colorbar;
xlabel('Productivity [kg {CO_2}/kg_s/day]', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('Electrical Energy [MJ/kg CO_2]', 'FontWeight', 'bold', 'FontSize', 20);
zlabel('g{H_2O}_{lost}/g{CO_2}', 'FontWeight', 'bold', 'FontSize', 20);
c = colorbar; c.Label.String = 'Electrical Energy [MJ/kg CO_2]'; c.Label.FontSize = 15; c.Label.FontWeight = 'bold'; 
set(gca, 'FontSize', 15, 'FontWeight', 'bold'); % Make axis markers bold
grid on; set(gca, 'LineWidth', 2); hold off;
% ylim([0 6]); xlim([0.05 0.25]); zlim([2 7])


%%

figure(3);
t = tiledlayout(2,4); % 2x4 grid
t.TileSpacing = 'compact'; % tighten spacing
t.Padding = 'compact';
% title(t,'Effects of Decision Variables on Energy and Productivity','FontWeight','bold','FontSize',14);

% --- Store axes handles ---
ax = gobjects(1,8);

% === Energy row (top, no x-labels/ticks) ===
ax(1) = nexttile(1);
scatter(v_ads, energy, 80, 'filled', 'MarkerFaceColor', [0.4940 0.1840 0.5560]);
ylabel('Electrical Energy [MJ/kg CO_2]','FontWeight','bold','FontSize',10);
set(gca,'XTickLabel',[]); % hide x ticks/labels

ax(2) = nexttile(2);
scatter(v_vs, energy, 80, 'filled', 'MarkerFaceColor', [0.8500 0.3250 0.0980]);
ylabel(''); set(gca,'YTickLabel',[]); 
set(gca,'XTickLabel',[]);

ax(3) = nexttile(3);
scatter(ads_time, energy, 80, 'filled', 'MarkerFaceColor', [0.4660 0.6740 0.1880]);
ylabel(''); set(gca,'YTickLabel',[]); 
set(gca,'XTickLabel',[]);

ax(4) = nexttile(4);
scatter(vs_time, energy, 80, 'filled', 'MarkerFaceColor', [0.3010 0.7450 0.9330]);
ylabel(''); set(gca,'YTickLabel',[]); 
set(gca,'XTickLabel',[]);

% === Productivity row (bottom, keep x-labels/ticks) ===
ax(5) = nexttile(5);
scatter(v_ads, prodtvy, 80, 'filled', 'MarkerFaceColor', [0.4940 0.1840 0.5560]);
xlabel('V_{ADS} [m/s]','FontWeight','bold','FontSize',10);
ylabel('Productivity [kg {CO_2}/kg_s/day]','FontWeight','bold','FontSize',10);

ax(6) = nexttile(6);
scatter(v_vs, prodtvy, 80, 'filled', 'MarkerFaceColor', [0.8500 0.3250 0.0980]);
xlabel('V_{VS} [m/s]','FontWeight','bold','FontSize',10);
ylabel(''); set(gca,'YTickLabel',[]);

ax(7) = nexttile(7);
scatter(ads_time, prodtvy, 80, 'filled', 'MarkerFaceColor', [0.4660 0.6740 0.1880]);
xlabel('t_{AD} [min]','FontWeight','bold','FontSize',10);
ylabel(''); set(gca,'YTickLabel',[]);

ax(8) = nexttile(8);
scatter(vs_time, prodtvy, 80, 'filled', 'MarkerFaceColor', [0.3010 0.7450 0.9330]);
xlabel('t_{VS} [min]','FontWeight','bold','FontSize',10);
ylabel(''); set(gca,'YTickLabel',[]);

% --- Link axes by row ---
linkaxes(ax(1:4),'y');  % Energy y-axis
linkaxes(ax(5:8),'y');  % Productivity y-axis

% --- Link axes by column ---
% linkaxes(ax(1:2:end),'x'); % v_ads (1,5), ads_time (3,7)
% linkaxes(ax(2:2:end),'x'); % v_vs (2,6), vs_time (4,8)

%%
figure(4); % Effects of variables on water loss
t = tiledlayout(2,4); % 2x4 grid
t.TileSpacing = 'compact'; % tighten spacing
t.Padding = 'compact';
% title(t,'Effects of Variables on Water Loss','FontWeight','bold','FontSize',14);

% --- Store axes handles ---
ax = gobjects(1,8);

% === Water Loss row (top, no x-labels/ticks) ===
ax(1) = nexttile(1);
scatter(v_ads, gH2Olost_gCO2, 80, 'filled', 'MarkerFaceColor', [0.4940 0.1840 0.5560]); % Purple
xlabel('V_{ADS} [m/s]','fontweight','bold','fontsize',10);
ylabel('g{H_2O} _{lost}/g{CO_2}','fontweight','bold','fontsize',10);
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);
set(gca,'XTickLabel',[]); % hide x ticks/labels

ax(2) = nexttile(2);
scatter(v_vs, gH2Olost_gCO2, 80, 'filled', 'MarkerFaceColor', [0.8500 0.3250 0.0980]); % Red-Orange
xlabel('V_{VS} [m/s]','fontweight','bold','fontsize',10);
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);
set(gca,'XTickLabel',[]); % hide x ticks/labels

ax(3) = nexttile(3);
scatter(ads_time, gH2Olost_gCO2, 80, 'filled', 'MarkerFaceColor', [0.4660 0.6740 0.1880]); % Green
xlabel('t_{AD} [min]','fontweight','bold','fontsize',10);
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);
set(gca,'XTickLabel',[]); % hide x ticks/labels

ax(4) = nexttile(4);
scatter(vs_time, gH2Olost_gCO2, 80, 'filled', 'MarkerFaceColor', [0.3010 0.7450 0.9330]); % Light Blue
xlabel('t_{VS} [min]','fontweight','bold','fontsize',10);
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);
set(gca,'XTickLabel',[]); % hide x ticks/labels

% === Water Used row (bottom, keep x-labels/ticks) ===
ax(5) = nexttile(5);
scatter(v_ads, gH2Oused_gCO2cap, 80, 'filled', 'MarkerFaceColor', [0.4940 0.1840 0.5560]); % Purple
xlabel('V_{ADS} [m/s]','fontweight','bold','fontsize',10);
ylabel('g{H_2O} _{proc}/g{CO_2}','fontweight','bold','fontsize',10);
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);

ax(6) = nexttile(6);
scatter(v_vs, gH2Oused_gCO2cap, 80, 'filled', 'MarkerFaceColor', [0.8500 0.3250 0.0980]); % Red-Orange
xlabel('V_{VS} [m/s]','fontweight','bold','fontsize',10);
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);

ax(7) = nexttile(7);
scatter(ads_time, gH2Oused_gCO2cap, 80, 'filled', 'MarkerFaceColor', [0.4660 0.6740 0.1880]); % Green
xlabel('t_{AD} [min]','fontweight','bold','fontsize',10);
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);

ax(8) = nexttile(8);
scatter(vs_time, gH2Oused_gCO2cap, 80, 'filled', 'MarkerFaceColor', [0.3010 0.7450 0.9330]); % Light Blue
xlabel('t_{VS} [min]','fontweight','bold','fontsize',10);
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);

% --- Link axes by row ---
linkaxes(ax(1:4),'y');  % Water loss y-axis
linkaxes(ax(5:8),'y');  % Water used y-axis

%% 
figure(7); % Times vs Velocities
subplot(2,1,1)
scatter(v_ads, ads_time, 80, 'filled', 'MarkerFaceColor', 'k'); % Black
xlabel('V_{AD} [m/s]','fontweight','bold','fontsize',10);
ylabel('t_{AD} [min]','fontweight','bold','fontsize',10);
hold on;
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);
% xlim([0 1]); ylim([40 160])
set(gca, 'YDir', 'normal'); % Set the y-axis direction to normal

% Add trend line
p = polyfit(v_ads, ads_time, 1); % Linear fit
yfit = polyval(p, v_ads);
plot(v_ads, yfit, 'r--', 'LineWidth', 2); % Trend line in red

subplot(2,1,2)
scatter(v_vs, vs_time, 80, 'filled', 'MarkerFaceColor', 'k'); % Black
xlabel('V_{VS} [m/s]','fontweight','bold','fontsize',10);
ylabel('t_{VS} [min]','fontweight','bold','fontsize',10);
% xlim([0 10]); ylim([40 180])
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);

% Add trend line
p = polyfit(v_vs, vs_time, 1); % Linear fit
yfit = polyval(p, v_vs);
hold on; % Ensure the trend line is added to the existing scatter plot
plot(v_vs, yfit, 'r--', 'LineWidth', 2); % Trend line in red

%% 
figure(8); % Water processed vs water lost
scatter(gH2Oused_gCO2cap, gH2Olost_gCO2, 80, 'filled', 'MarkerFaceColor', 'k'); % Black
xlabel('g{H_2O} _{proc}/g{CO_2}','fontweight','bold','fontsize',10);
ylabel('g{H_2O} _{lost}/g{CO_2}','fontweight','bold','fontsize',10);
hold on;
% ylim([2 7]); xlim([70 3500]) % Updated limits
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);

%% Bubble Charts
% Calculate energy range
energy_min = min(energy);
energy_max = max(energy);

prodtvy_min = min(prodtvy);
prodtvy_max = max(prodtvy);

% Create bubble chart
figure(9)
bubblechart(v_vs, v_ads, ads_time, energy, 'MarkerFaceAlpha', 0.7)
xlabel('v_{VS} [m/s]','fontweight','bold','fontsize',12); 
ylabel('v_{AD} [m/s]','fontweight','bold','fontsize',12)
colorbar; colormap(flipud(autumn)) % or use 'turbo', 'hot', etc.
c = colorbar; set(c.Label, 'String', 'Energy [MJ/kg CO_2]', 'FontSize', 12, 'FontWeight', 'bold'); 
set(gca, 'LineWidth', 2); % Make plot boundaries bold
time_range_str = sprintf('Time Range: %.2f–%.2f min', min(ads_time), max(ads_time));
text(min(v_vs) + 0.6*range(v_vs), min(v_ads) , time_range_str, ...
    'fontweight','bold','fontsize',12);
 % xlim([-0.5 10.5]); ylim([-0.5 1.1])


figure(10)
bubblechart(v_vs, v_ads, vs_time, energy, ...
    'MarkerFaceAlpha', 0.7)
xlabel('v_{VS} [m/s]','fontweight','bold','fontsize',12); 
ylabel('v_{AD} [m/s]','fontweight','bold','fontsize',12)
colorbar; colormap(flipud(autumn)) % or use 'turbo', 'hot', etc.
c = colorbar; set(c.Label, 'String', 'Energy [MJ/kg CO_2]', 'FontSize', 12, 'FontWeight', 'bold'); 
set(gca, 'LineWidth', 2); % Make plot boundaries bold
vs_time_range_str = sprintf('Time Range: %.2f–%.2f min', min(vs_time), max(vs_time));
text(min(v_vs) + 0.6*range(v_vs), min(v_ads) + 0.02*range(v_ads), vs_time_range_str, ...
    'fontweight','bold','fontsize',12);
 % xlim([-0.5 10.5]); ylim([-0.5 1.1])


figure(11)
bubblechart(v_vs, v_ads, ads_time, prodtvy, ...
    'MarkerFaceAlpha', 0.7)
xlabel('v_{VS} [m/s]','fontweight','bold','fontsize',12); 
ylabel('v_{AD} [m/s]','fontweight','bold','fontsize',12)
colorbar; colormap(flipud(summer)) % or use 'turbo', 'hot', etc.
c = colorbar; set(c.Label, 'String', 'Productivity [kg {CO_2}/kg_s/day]', 'FontSize', 12, 'FontWeight', 'bold'); 
set(gca, 'LineWidth', 2); % Make plot boundaries bold
% ylim([65 85]); xlim([115 125]);
ads_time_range_str = sprintf('Adsorption Time Range: %.2f–%.2f min', min(ads_time), max(ads_time));
text(min(v_vs), max(v_ads) + 0.02*range(v_ads), ads_time_range_str, ...
    'fontweight','bold','fontsize',12);

figure(12)
bubblechart(vs_time, ads_time, vs_time, v_vs, ...
    'MarkerFaceAlpha', 0.7)
xlabel('t_{VS} [min]','fontweight','bold','fontsize',12); 
ylabel('t_{AD} [min]','fontweight','bold','fontsize',12)
colorbar; colormap(flipud(summer)) % or use 'turbo', 'hot', etc.
c = colorbar; set(c.Label, 'String', 'v_{VS} [m/s]', 'FontSize', 12, 'FontWeight', 'bold'); 
set(gca, 'LineWidth', 2); % Make plot boundaries bold
% ylim([65 85]); xlim([115 125]);
vs_time_range_str = sprintf('Time Range: %.2f–%.2f min', min(vs_time), max(vs_time));
text(min(vs_time), max(ads_time) + 0.02*range(ads_time), vs_time_range_str, ...
    'fontweight','bold','fontsize',12);


figure(13)
bubblechart(v_ads, energy, gH2Olost_gCO2);%, ...
    % 'MarkerFaceAlpha', 0.7)
xlabel('v_{AD} [m/s]','fontweight','bold','fontsize',12); 
ylabel('Electrical Energy [MJ/kg CO_2]','fontweight','bold','fontsize',12)
colorbar; colormap(sky) 
c = colorbar; set(c.Label, 'String', 'g_{H_2O_{lost}} / g_{CO_2}', 'FontSize', 12, 'FontWeight', 'bold'); 
set(gca, 'LineWidth', 2); % Make plot boundaries bold
% ylim([0.5 2.5]); xlim([0 1]);


figure(14)
bubblechart(v_vs,prodtvy, gH2Olost_gCO2, ...
    'MarkerFaceAlpha', 0.7)
xlabel('v_{VS} [m/s]','fontweight','bold','fontsize',12); 
ylabel('Productivity [kg {CO_2}/kg_s/day]','fontweight','bold','fontsize',12)
colorbar; colormap(sky) 
c = colorbar; set(c.Label, 'String', 'g_{H_2O_{lost}} / g_{CO_2}', 'FontSize', 12, 'FontWeight', 'bold'); 
set(gca, 'LineWidth', 2); % Make plot boundaries bold
% ylim([0.5 2.5]); xlim([0.9 1.02]);



%% Bubble Charts
% % Calculate energy range
% energy_min = min(energy);
% energy_max = max(energy);
% 
% prodtvy_min = min(prodtvy);
% prodtvy_max = max(prodtvy);
% 
% % Create bubble chart
% figure(5)
% bubblechart(vs_time, ads_time, energy, v_ads,'MarkerFaceAlpha', 0.7)
% xlabel('t_{VS} [min]','fontweight','bold','fontsize',12); 
% ylabel('t_{AD} [min]','fontweight','bold','fontsize',12)
% colorbar; colormap(parula) % or use 'turbo', 'hot', etc.
% c = colorbar; set(c.Label, 'String', 'v_{AD} [m/s]', 'FontSize', 12, 'FontWeight', 'bold'); 
% % ylim([65 85]); xlim([115 125]);
% set(gca, 'LineWidth', 2); % Make plot boundaries bold
% energy_range_str = sprintf('Energy Range: %.2f–%.2f MJ/kg CO_2', energy_min, energy_max);
% text(min(vs_time), max(ads_time) + 0.02*range(ads_time), energy_range_str, ...
%     'fontweight','bold','fontsize',12);
% xlim([60 180]); ylim([40 160])
% 
% 
% 
% figure(6)
% bubblechart(vs_time, ads_time, energy, v_vs, ...
%     'MarkerFaceAlpha', 0.7)
% xlabel('t_{VS} [min]','fontweight','bold','fontsize',12); 
% ylabel('t_{AD} [min]','fontweight','bold','fontsize',12)
% colorbar; colormap(parula) % or use 'turbo', 'hot', etc.
% c = colorbar; set(c.Label, 'String', 'v_{VS} [m/s]', 'FontSize', 12, 'FontWeight', 'bold'); 
% % ylim([65 85]); xlim([115 125]);
% set(gca, 'LineWidth', 2); % Make plot boundaries bold
% energy_range_str = sprintf('Energy Range: %.2f–%.2f MJ/kg CO_2', energy_min, energy_max);
% text(min(vs_time), max(ads_time) + 0.02*range(ads_time), energy_range_str, ...
%     'fontweight','bold','fontsize',12);
% xlim([60 180]); ylim([40 160])
% 
% 
% 
% figure(7)
% bubblechart(vs_time, ads_time, prodtvy, v_ads, ...
%     'MarkerFaceAlpha', 0.7)
% xlabel('t_{VS} [min]','fontweight','bold','fontsize',12); 
% ylabel('t_{AD} [min]','fontweight','bold','fontsize',12)
% colorbar; colormap(parula) % or use 'turbo', 'hot', etc.
% c = colorbar; set(c.Label, 'String', 'v_{AD} [m/s]', 'FontSize', 12, 'FontWeight', 'bold'); 
% set(gca, 'LineWidth', 2); % Make plot boundaries bold
% % ylim([65 85]); xlim([115 125]);
% prodtvy_range_str = sprintf('Productivity Range: %.2f–%.2f [kg {CO_2}/kg_s/day]', prodtvy_min, prodtvy_max);
% text(min(vs_time), max(ads_time) + 0.02*range(ads_time), prodtvy_range_str, ...
%     'fontweight','bold','fontsize',12);
% xlim([60 180]); ylim([40 160])
% 
% 
% figure(8)
% bubblechart(vs_time, ads_time, prodtvy, v_vs, ...
%     'MarkerFaceAlpha', 0.7)
% xlabel('t_{VS} [min]','fontweight','bold','fontsize',12); 
% ylabel('t_{AD} [min]','fontweight','bold','fontsize',12)
% colorbar; colormap(parula) % or use 'turbo', 'hot', etc.
% c = colorbar; set(c.Label, 'String', 'v_{VS} [m/s]', 'FontSize', 12, 'FontWeight', 'bold'); 
% set(gca, 'LineWidth', 2); % Make plot boundaries bold
% % ylim([65 85]); xlim([115 125]);
% prodtvy_range_str = sprintf('Productivity Range: %.2f–%.2f [kg {CO_2}/kg_s/day]', prodtvy_min, prodtvy_max);
% text(min(vs_time), max(ads_time) + 0.02*range(ads_time), prodtvy_range_str, ...
%     'fontweight','bold','fontsize',12);
% xlim([60 180]); ylim([40 160])
%%
figure(20);
t = tiledlayout(3,4); % 3x4 grid
t.TileSpacing = 'compact'; % Increase space between tiles
t.Padding = 'compact';
% title(t,'Effects of Decision Variables on Energy, Productivity, and Water Lost','FontWeight','bold','FontSize',14);

ax = gobjects(1,12); % 12 subplots total
% === Row 1: Energy ===
ax(1) = nexttile(1);
scatter(v_ads, energy, 80, 'filled', 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
ylabel('E.Energy [MJ/kg CO_2]','FontWeight','bold','FontSize',12);
set(gca,'XTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold
% ylim([0 10]); % Set y limits

ax(2) = nexttile(2);
scatter(v_vs, energy, 80, 'filled', 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
ylabel(''); set(gca,'YTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca,'XTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold
% ylim([0 10]); % Set y limits

ax(3) = nexttile(3);
scatter(ads_time, energy, 80, 'filled', 'MarkerFaceColor',[0.4660 0.6740 0.1880]);
ylabel(''); set(gca,'YTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca,'XTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold
% ylim([0 10]); % Set y limits

ax(4) = nexttile(4);
scatter(vs_time, energy, 80, 'filled', 'MarkerFaceColor',[0.3010 0.7450 0.9330]);
ylabel(''); set(gca,'YTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca,'XTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold
% ylim([0 10]); % Set y limits

% === Row 2: Productivity ===
ax(5) = nexttile(5);
scatter(v_ads, prodtvy, 80, 'filled', 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
ylabel('Prod [kg {CO_2}/kg_s/day]','FontWeight','bold','FontSize',12);
set(gca,'XTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold
% ylim([0.05 0.15]); % Set y limits

ax(6) = nexttile(6);
scatter(v_vs, prodtvy, 80, 'filled', 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
ylabel(''); set(gca,'YTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca,'XTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold
% ylim([0.05 0.15]); % Set y limits

ax(7) = nexttile(7);
scatter(ads_time, prodtvy, 80, 'filled', 'MarkerFaceColor',[0.4660 0.6740 0.1880]);
ylabel(''); set(gca,'YTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca,'XTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold
% ylim([0.05 0.15]); % Set y limits

ax(8) = nexttile(8);
scatter(vs_time, prodtvy, 80, 'filled', 'MarkerFaceColor',[0.3010 0.7450 0.9330]);
ylabel(''); set(gca,'YTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca,'XTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold
% ylim([0.05 0.15]); % Set y limits

% === Row 3: Water Lost (keep x-axis labels here) ===
ax(9) = nexttile(9);
scatter(v_ads, gH2Olost_gCO2, 80, 'filled', 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
xlabel('V_{ADS} [m/s]','FontWeight','bold','FontSize',12);
ylabel('g_{H_2O}/g_{CO_2} lost','FontWeight','bold','FontSize',12);
set(gca, 'LineWidth', 2, 'FontWeight', 'bold', 'FontSize', 12); % Make axis marker bold

ax(10) = nexttile(10);
scatter(v_vs, gH2Olost_gCO2, 80, 'filled', 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
xlabel('V_{VS} [m/s]','FontWeight','bold','FontSize',12);
ylabel(''); set(gca,'YTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold

ax(11) = nexttile(11);
scatter(ads_time, gH2Olost_gCO2, 80, 'filled', 'MarkerFaceColor',[0.4660 0.6740 0.1880]);
xlabel('t_{AD} [min]','FontWeight','bold','FontSize',12);
ylabel(''); set(gca,'YTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold

ax(12) = nexttile(12);
scatter(vs_time, gH2Olost_gCO2, 80, 'filled', 'MarkerFaceColor',[0.3010 0.7450 0.9330]);
xlabel('t_{VS} [min]','FontWeight','bold','FontSize',12);
ylabel(''); set(gca,'YTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold

% --- Link axes ---
linkaxes(ax(1:4),'y');   % Energy row
linkaxes(ax(5:8),'y');   % Productivity row
linkaxes(ax(9:12),'y');  % Water Lost row

% Link x by column (align decision variables)
linkaxes([ax(1) ax(5) ax(9)],'x'); % v_ads
linkaxes([ax(2) ax(6) ax(10)],'x'); % v_vs
linkaxes([ax(3) ax(7) ax(11)],'x'); % ads_time
linkaxes([ax(4) ax(8) ax(12)],'x'); % vs_time

% Set x limits for ads_time and vs_time
% xlim(ax(3), [0 120]); % Set x limits for ads_time
% xlim(ax(2), [5 10]); % Set x limits for vs_velocity
