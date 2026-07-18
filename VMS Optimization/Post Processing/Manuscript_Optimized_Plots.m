% Define blue color
blueColor = [0 0.4470 0.7410];

%% Figure 1: Energy vs Productivity
% Convert units for manuscript
energy   = energy * (3.6/1000)           ; % from kWh/tonCO2 to MJ/kgCO2
prodtvy  = prodtvy *((24*0.044)/(1060))  ; % from mol/m3/hr to kgCO2/kgs/day
cs_time  = cs_time/60                    ; % convert sec to min
vs_time  = vs_time/60                    ; % convert sec to min
des_time = des_time/60                   ; % convert sec to min

% Combine all variables into a table
T = table(energy(:), prodtvy(:), cs_time(:), vs_time(:), des_time(:), ...
          vair(:), vvs(:), gH2Olost_gCO2(:), gH2Oused_gCO2cap(:), CO2_Purity(:), ...
          'VariableNames', {'Energy', 'Productivity', 'AdsTime', 'VSTime', 'DSTime', ...
                            'vair', 'vvs', 'gH2OLost_gCO2', 'gH2OUsed_gCO2Cap', 'CO2 Purity'});
writetable(T, 'Optimized_3obj_P50G30_50RH_MK_EnCorrection_fullyloaded.xlsx');

figure(1); clf;
scatter3(prodtvy, energy, gH2Olost_gCO2, 150, gH2Olost_gCO2, 'filled'); colormap(sky); colorbar;
xlabel('Productivity [kg {CO_2}/kg_s/day]', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('Electrical Energy [MJ/kg CO_2]', 'FontWeight', 'bold', 'FontSize', 20);
zlabel('g{H_2O}_{lost}/g{CO_2}', 'FontWeight', 'bold', 'FontSize', 20);
c = colorbar; c.Label.String = 'g{H_2O}_{lost}/g{CO_2}'; c.Label.FontSize = 15; c.Label.FontWeight = 'bold'; 
set(gca, 'FontSize', 15, 'FontWeight', 'bold'); % Make axis markers bold
grid on; set(gca, 'LineWidth', 2); hold off;
% ylim([0 6]); 
% xlim([0.0 1]); 
% zlim([2 7])

figure(2); clf;
scatter3(prodtvy, energy, gH2Olost_gCO2, 150, energy, 'filled'); colormap(flipud(autumn)); colorbar;
xlabel('Productivity [kg {CO_2}/kg_s/day]', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('Electrical Energy [MJ/kg CO_2]', 'FontWeight', 'bold', 'FontSize', 20);
zlabel('g{H_2O}_{lost}/g{CO_2}', 'FontWeight', 'bold', 'FontSize', 20);
c = colorbar; c.Label.String = 'Electrical Energy [MJ/kg CO_2]'; c.Label.FontSize = 15; c.Label.FontWeight = 'bold'; 
set(gca, 'FontSize', 15, 'FontWeight', 'bold'); % Make axis markers bold
grid on; set(gca, 'LineWidth', 2); hold off;
% ylim([0 6]); xlim([0.05 0.25]); zlim([2 7])
% xlim([0.0 1]); 


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
scatter(vair, energy, 80, 'filled', 'MarkerFaceColor', [0.4940 0.1840 0.5560]);
ylabel('Electrical Energy [MJ/kg CO_2]','FontWeight','bold','FontSize',10);
set(gca,'XTickLabel',[]); % hide x ticks/labels

ax(2) = nexttile(2);
scatter(vvs, energy, 80, 'filled', 'MarkerFaceColor', [0.8500 0.3250 0.0980]);
ylabel(''); set(gca,'YTickLabel',[]); 
set(gca,'XTickLabel',[]);

ax(3) = nexttile(3);
scatter(cs_time, energy, 80, 'filled', 'MarkerFaceColor', [0.4660 0.6740 0.1880]);
ylabel(''); set(gca,'YTickLabel',[]); 
set(gca,'XTickLabel',[]);

ax(4) = nexttile(4);
scatter(vs_time, energy, 80, 'filled', 'MarkerFaceColor', [0.3010 0.7450 0.9330]);
ylabel(''); set(gca,'YTickLabel',[]); 
set(gca,'XTickLabel',[]);

% === Productivity row (bottom, keep x-labels/ticks) ===
ax(5) = nexttile(5);
scatter(vair, prodtvy, 80, 'filled', 'MarkerFaceColor', [0.4940 0.1840 0.5560]);
xlabel('V_{ADS} [m/s]','FontWeight','bold','FontSize',10);
ylabel('Productivity [kg {CO_2}/kg_s/day]','FontWeight','bold','FontSize',10);

ax(6) = nexttile(6);
scatter(vvs, prodtvy, 80, 'filled', 'MarkerFaceColor', [0.8500 0.3250 0.0980]);
xlabel('V_{VS} [m/s]','FontWeight','bold','FontSize',10);
ylabel(''); set(gca,'YTickLabel',[]);

ax(7) = nexttile(7);
scatter(cs_time, prodtvy, 80, 'filled', 'MarkerFaceColor', [0.4660 0.6740 0.1880]);
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
% linkaxes(ax(1:2:end),'x'); % vair (1,5), ads_time (3,7)
% linkaxes(ax(2:2:end),'x'); % vvs (2,6), vs_time (4,8)

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
scatter(vair, gH2Olost_gCO2, 80, 'filled', 'MarkerFaceColor', [0.4940 0.1840 0.5560]); % Purple
xlabel('V_{AIR} [m/s]','fontweight','bold','fontsize',10);
ylabel('g{H_2O} _{lost}/g{CO_2}','fontweight','bold','fontsize',10);
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);
set(gca,'XTickLabel',[]); % hide x ticks/labels

ax(2) = nexttile(2);
scatter(vvs, gH2Olost_gCO2, 80, 'filled', 'MarkerFaceColor', [0.8500 0.3250 0.0980]); % Red-Orange
xlabel('V_{VS} [m/s]','fontweight','bold','fontsize',10);
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);
set(gca,'XTickLabel',[]); % hide x ticks/labels

ax(3) = nexttile(3);
scatter(cs_time, gH2Olost_gCO2, 80, 'filled', 'MarkerFaceColor', [0.4660 0.6740 0.1880]); % Green
xlabel('t_{CS} [min]','fontweight','bold','fontsize',10);
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);
set(gca,'XTickLabel',[]); % hide x ticks/labels

ax(4) = nexttile(4);
scatter(vs_time, gH2Olost_gCO2, 80, 'filled', 'MarkerFaceColor', [0.3010 0.7450 0.9330]); % Light Blue
xlabel('t_{VS} [min]','fontweight','bold','fontsize',10);
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);
set(gca,'XTickLabel',[]); % hide x ticks/labels

% === Water Used row (bottom, keep x-labels/ticks) ===
ax(5) = nexttile(5);
scatter(vair, gH2Oused_gCO2cap, 80, 'filled', 'MarkerFaceColor', [0.4940 0.1840 0.5560]); % Purple
xlabel('V_{AIR} [m/s]','fontweight','bold','fontsize',10);
ylabel('g{H_2O} _{proc}/g{CO_2}','fontweight','bold','fontsize',10);
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);

ax(6) = nexttile(6);
scatter(vvs, gH2Oused_gCO2cap, 80, 'filled', 'MarkerFaceColor', [0.8500 0.3250 0.0980]); % Red-Orange
xlabel('V_{VS} [m/s]','fontweight','bold','fontsize',10);
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);

ax(7) = nexttile(7);
scatter(cs_time, gH2Oused_gCO2cap, 80, 'filled', 'MarkerFaceColor', [0.4660 0.6740 0.1880]); % Green
xlabel('t_{CS} [min]','fontweight','bold','fontsize',10);
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);

ax(8) = nexttile(8);
scatter(vs_time, gH2Oused_gCO2cap, 80, 'filled', 'MarkerFaceColor', [0.3010 0.7450 0.9330]); % Light Blue
xlabel('t_{VS} [min]','fontweight','bold','fontsize',10);
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);

% --- Link axes by row ---
linkaxes(ax(1:4),'y');  % Water loss y-axis
linkaxes(ax(5:8),'y');  % Water used y-axis

%% 
figure(5); % Times vs Velocities
subplot(2,1,1)
scatter(vair, cs_time, 80, 'filled', 'MarkerFaceColor', 'k'); % Black
xlabel('V_{AIR} [m/s]','fontweight','bold','fontsize',10);
ylabel('t_{CS} [min]','fontweight','bold','fontsize',10);
hold on;
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);
% xlim([0 1]); ylim([40 160])
set(gca, 'YDir', 'normal'); % Set the y-axis direction to normal

% Add trend line
% p = polyfit(vair, cs_time, 1); % Linear fit
% yfit = polyval(p, vair);
% plot(vair, yfit, 'r--', 'LineWidth', 2); % Trend line in red

subplot(2,1,2)
scatter(vvs, vs_time, 80, 'filled', 'MarkerFaceColor', 'k'); % Black
xlabel('V_{VS} [m/s]','fontweight','bold','fontsize',10);
ylabel('t_{VS} [min]','fontweight','bold','fontsize',10);
% xlim([0 10]); ylim([40 180])
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); set(gca, 'LineWidth', 2);



%%
figure(6); %Effects of Decision Variables on Energy, Productivity, and Water Lost
t = tiledlayout(3,4); % 3x4 grid
t.TileSpacing = 'compact'; % Increase space between tiles
t.Padding = 'compact';

idx_vvs = vvs<5;
ax = gobjects(1,12); % 12 subplots total
% === Row 1: Energy ===
ax(1) = nexttile(1);
scatter(vair(idx_vvs), energy(idx_vvs), 80, 'filled', 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
ylabel({'Electrical Energy'; '[MJ/kg CO_2]'}, 'FontWeight', 'bold', 'FontSize', 12);
set(gca,'XTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold
ylim([0 30]); % Set y limits

ax(2) = nexttile(2);
scatter(vvs(idx_vvs), energy(idx_vvs), 80, 'filled', 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
ylabel(''); set(gca,'YTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca,'XTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold
% ylim([0 20]); % Set y limits

ax(3) = nexttile(3);
scatter(cs_time(idx_vvs), energy(idx_vvs), 80, 'filled', 'MarkerFaceColor',[0.4660 0.6740 0.1880]);
ylabel(''); set(gca,'YTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca,'XTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold
% ylim([0 10]); % Set y limits

ax(4) = nexttile(4);
scatter(vs_time(idx_vvs), energy(idx_vvs), 80, 'filled', 'MarkerFaceColor',[0.3010 0.7450 0.9330]);
ylabel(''); set(gca,'YTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca,'XTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold
% ylim([0 10]); % Set y limits

% === Row 2: Productivity ===
ax(5) = nexttile(5);
scatter(vair(idx_vvs), prodtvy(idx_vvs), 80, 'filled', 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
ylabel({'Productivity'; '[kg {CO_2}/kg_s/day]'},'FontWeight','bold','FontSize',12);
set(gca,'XTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold
ylim([0 1]); % Set y limits

ax(6) = nexttile(6);
scatter(vvs(idx_vvs), prodtvy(idx_vvs), 80, 'filled', 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
ylabel(''); set(gca,'YTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca,'XTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold
% ylim([0.05 0.15]); % Set y limits

ax(7) = nexttile(7);
scatter(cs_time(idx_vvs), prodtvy(idx_vvs), 80, 'filled', 'MarkerFaceColor',[0.4660 0.6740 0.1880]);
ylabel(''); set(gca,'YTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca,'XTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold
% ylim([0.05 0.15]); % Set y limits

ax(8) = nexttile(8);
scatter(vs_time(idx_vvs), prodtvy(idx_vvs), 80, 'filled', 'MarkerFaceColor',[0.3010 0.7450 0.9330]);
ylabel(''); set(gca,'YTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca,'XTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold
% ylim([0.05 0.15]); % Set y limits

% === Row 3: Water Lost (keep x-axis labels here) ===
ax(9) = nexttile(9);
scatter(vair(idx_vvs), gH2Olost_gCO2(idx_vvs), 80, 'filled', 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
xlabel('V_{AIR} [m/s]','FontWeight','bold','FontSize',12);
ylabel('g_{H_2O}/g_{CO_2} lost','FontWeight','bold','FontSize',12);
set(gca, 'LineWidth', 2, 'FontWeight', 'bold', 'FontSize', 12); % Make axis marker bold
ylim([0 5]);

ax(10) = nexttile(10);
scatter(vvs(idx_vvs), gH2Olost_gCO2(idx_vvs), 80, 'filled', 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
xlabel('V_{VS} [m/s]','FontWeight','bold','FontSize',12);
ylabel(''); set(gca,'YTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold

ax(11) = nexttile(11);
scatter(cs_time(idx_vvs), gH2Olost_gCO2(idx_vvs), 80, 'filled', 'MarkerFaceColor',[0.4660 0.6740 0.1880]);
xlabel('t_{CS} [min]','FontWeight','bold','FontSize',12);
ylabel(''); set(gca,'YTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold

ax(12) = nexttile(12);
scatter(vs_time(idx_vvs), gH2Olost_gCO2(idx_vvs), 80, 'filled', 'MarkerFaceColor',[0.3010 0.7450 0.9330]);
xlabel('t_{VS} [min]','FontWeight','bold','FontSize',12);
ylabel(''); set(gca,'YTickLabel',[], 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'LineWidth', 2); % Make axis marker bold

% --- Link axes ---
linkaxes(ax(1:4),'y');   % Energy row
linkaxes(ax(5:8),'y');   % Productivity row
linkaxes(ax(9:12),'y');  % Water Lost row

% Link x by column (align decision variables)
linkaxes([ax(1) ax(5) ax(9)],'x'); % vair
linkaxes([ax(2) ax(6) ax(10)],'x'); % vvs
linkaxes([ax(3) ax(7) ax(11)],'x'); % cs_time
linkaxes([ax(4) ax(8) ax(12)],'x'); % vs_time

% Set x limits for vair, vvs, ads_time, and vs_time
xlim(ax(5), [0 1]);    % Set x limits for vair
xlim(ax(6), [0 5]);   % Set x limits for vvs
xlim(ax(3), [0 200]); % Set x limits for cs_time
xlim(ax(4), [0 200]); % Set x limits for vs_time
