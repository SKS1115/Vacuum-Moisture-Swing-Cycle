

function PerformanceMetrics = MechanicalEnergy(a, b, c, d, e, t1, t2, t3, t4, t5, RH, N, L, ColumnParam, CondParam, Times)

    r_in    = ColumnParam(3) ;
    eb      = ColumnParam(5) ;      
    rho_s   = CondParam(1)   ;    
    gamma   = CondParam(9)   ;
    R       = CondParam(14)  ;   
    eta     = CondParam(15)  ;     
    PH      = CondParam(20)  ;    
    Tsat    = CondParam(24)  ;

    Psat_vs = PressureVaporSat(Tsat) ;               

% Calculate Sorbent Properties
    CycleTime   = sum(Times)                 ;   % Total time of cycle
    Vol_Sorbent = (1-eb)*pi*(r_in)^2*L       ;   % Volume of sorbent (m3)
    MM_CO2      = 44.01                         ;   % Molar mass of CO2 (g/mol)
    MM_H2O      = 18.02                         ;   % Molar mass of H2O (g/mol)
    g_sorbent   = (rho_s*Vol_Sorbent)*1000   ;   % Mass of sorbent (grams)
    gam1        = gamma/ (gamma-1);
    gam2        = (gamma-1)/gamma;

% Extract Data Function
    function [yCO2, yH2O, Press, Temp, qCO2, qH2O, yN2] = extract_params(matrix,N)       
        yCO2 = matrix(:, 1:N+2);
        yH2O = matrix(:, N+3:2*N+4);
        Press= matrix(:, 2*N+5:3*N+6);
        Temp = matrix(:, 5*N+11:6*N+12);
        qCO2 = matrix(:, 3*N+7:4*N+8);
        qH2O = matrix(:, 4*N+9:5*N+10);
        yN2  = max(1 - (yCO2 + yH2O), 0); % Ensure no negative values
    end

    % Extract Data
    [yCO2_A, yH2O_A, Press_A, Temp_A, qCO2_A, qH2O_A, yN2_A] = extract_params(a,N);
    [yCO2_B, yH2O_B, Press_B, Temp_B, qCO2_B, qH2O_B, yN2_B] = extract_params(b,N);
    [yCO2_C, yH2O_C, Press_C, Temp_C, qCO2_C, qH2O_C, yN2_C] = extract_params(c,N);
    [yCO2_D, yH2O_D, Press_D, Temp_D, qCO2_D, qH2O_D, yN2_D] = extract_params(d,N);
    [yCO2_E, yH2O_E, Press_E, Temp_E, qCO2_E, qH2O_E, yN2_E] = extract_params(e,N);


%%%.................Extract Key Inputs from Mass Balance Function............................%%%
MassBalance     = MassBalanceFxn(a,b,c,d,e,t1,t2,t3,t4,t5,N,L,ColumnParam,CondParam);
CO2_molarflow   = MassBalance.CO2_molarflow ; 
CO2_moles       = MassBalance.CO2_moles     ;
H2O_molarflow   = MassBalance.H2O_molarflow ; 
H2O_moles       = MassBalance.H2O_moles     ; 
N2_molarflow    = MassBalance.N2_molarflow  ; 
N2_moles        = MassBalance.N2_moles      ; 

nCO2in_pz   = CO2_molarflow{1}; nH2Oin_pz   = H2O_molarflow{1}; nN2in_pz   = N2_molarflow{1}; 
nCO2in_cs   = CO2_molarflow{2}; nH2Oin_cs   = H2O_molarflow{2}; nN2in_cs   = N2_molarflow{2}; 
nCO2out_ev  = CO2_molarflow{3}; nH2Oout_ev  = H2O_molarflow{3}; nN2out_ev  = N2_molarflow{3}; 
nCO2out_vs  = CO2_molarflow{4}; nH2Oout_vs  = H2O_molarflow{4}; nN2out_vs  = N2_molarflow{4}; 
nCO2out_des = CO2_molarflow{5}; nH2Oout_des = H2O_molarflow{5}; nN2out_des = N2_molarflow{5}; 

molCO2in_pz  = CO2_moles{1}; molH2Oin_pz  = H2O_moles{1}; molH2Oout_des = H2O_moles{5};
molCO2in_cs  = CO2_moles{2}; molH2Oin_cs  = H2O_moles{2}; molN2out_vs   = N2_moles{1}; 
molCO2out_cs = CO2_moles{3}; molH2Oin_vs  = H2O_moles{3}; molN2out_des  = N2_moles{2}; 
molCO2out_vs = CO2_moles{4}; molH2Oout_cs = H2O_moles{4}; 
molCO2out_des= CO2_moles{5}; molH2Oout_vs = H2O_moles{5};  

% Flow compositions
    Total_molarIn_cs = trapz(t2, (nCO2in_cs + nH2Oin_cs + nN2in_cs)); % mol/s

% Calculate Key Performance metrics
    molCO2_sorbed    = molCO2in_pz  + (molCO2in_cs-molCO2out_cs)   ;   % moles of CO2 adsorbed
    gramCO2_sorbed   = molCO2_sorbed*MM_CO2                       ;
    molCO2_captured    = molCO2out_vs + molCO2out_des                ;   % moles of CO2 captured
    gramCO2_captured   = molCO2_captured*MM_CO2                       ;   % mass of CO2 captured (grams)
    tonCO2_captured    = gramCO2_captured/1e6                        ;   % mass of CO2 captured (tons)
    gCO2_gSorbent      = gramCO2_captured/g_sorbent                  ;   % grams CO2/grams Sorbent
    CO2purity          = ((molCO2_captured)/(molCO2out_vs + molN2out_vs + molCO2out_des + molN2out_des))*100 ;
    CO2Recovery        = (molCO2_captured/molCO2_sorbed)*100;

    molH2O_used        = molH2Oin_vs                     ;
    molH2O_recovered   = molH2Oout_vs + molH2Oout_des    ; 
    molH2O_lost        = molH2Oout_cs - (molH2Oin_cs + molH2Oin_pz) ;  
    H2ORecovery        = (molH2O_recovered/molH2O_used)*100    ; 
    
    gH2O_used      = molH2O_used*MM_H2O                 ;  % grams of H2O used 
    gH2O_recovered = molH2O_recovered*MM_H2O            ;  % grams of H2O recovered 
    gH2O_lost      = molH2O_lost*MM_H2O                 ;  % grams of H2O lost during adsorption
    tonsH2O_recovered  = gH2O_recovered/1e6            ;  % tons of H2O recovered 
    gH2Oused_gCO2   = gH2O_used/gramCO2_captured    ;  % grams H2O recovered/grams CO2 captured
    gH2Olost_gCO2   = gH2O_lost/gramCO2_captured    ;
 
%%%................... Total Energy Consumption...................

% Energy Consumption for Each Step.
    
% 1) Pressurization (Inlet)
    mflowin_pz    = (nCO2in_pz + nH2Oin_pz + nN2in_pz);
    Energydt_pz   = (1/eta)*(gam1).*(mflowin_pz.*R.*Temp_A(:,1)).*(((Press_A(:,1)./PH).^gam2)-1) ;
    EnergykWh_pz  = (trapz(t1, Energydt_pz))/3.6e6    ; % Blower energy in kWh

% 2) Blower Energy during Adsorption 
    mflowin_cs    = (nCO2in_cs + nH2Oin_cs + nN2in_cs);
    Energydt_cs   = ((1/eta)*(gam1).*(mflowin_cs.*R.*Temp_B(:,1)).*(((Press_B(:,1)./PH).^gam2)-1)) ; 
    EnergykWh_cs  = (trapz(t2, Energydt_cs))/3.6e6    ; % Blower energy in kWh
    EnergykW_cs = EnergykWh_cs/(t2(end)/3600);

% 3) Vacuum Energy during Air Evacuation
    mflowout_ev    = (nCO2out_ev + nH2Oout_ev + nN2out_ev);
    Energydt_ev    = (1/eta)*(gam1).*(mflowout_ev.*R.*Temp_C(:,end)).*(((PH./Press_C(:,end)).^gam2)-1) ;
    EnergykWh_ev   = (trapz(t3, Energydt_ev))/3.6e6    ; % Vacuum energy in kWh
    EnergykW_ev    = EnergykWh_ev/(t3(end)/3600);

% 4) Vacuum Energy during Vapor Stripping
    Edt_CO2_vs     = (1/eta)*(gam1).*(nCO2out_vs.*R.*Temp_D(:,end)).*(((PH./Press_D(:,end)).^gam2)-1) ;
    Edt_H2O_vs     = (1/eta)*(gam1).*(nH2Oout_vs.*R.*Temp_D(:,end)).*(((Psat_vs./Press_D(:,end)).^gam2)-1) ;
    Edt_N2_vs      = (1/eta)*(gam1).*(nN2out_vs.*R.*Temp_D(:,end)).*(((PH./Press_D(:,end)).^gam2)-1) ;
    EnergykWh_vs   = (trapz(t4, (Edt_CO2_vs + Edt_H2O_vs + Edt_N2_vs)))/3.6e6;

% 5)  Vacuum Energy during final desorption
    Edt_CO2_des    = (1/eta)*(gam1).*(nCO2out_des.*R.*Temp_E(:,end)).*(((PH./Press_E(:,end)).^gam2)-1) ;
    Edt_H2O_des    = (1/eta)*(gam1).*(nH2Oout_des.*R.*Temp_E(:,end)).*(((Psat_vs./Press_E(:,end)).^gam2)-1) ;
    Edt_N2_des     = (1/eta)*(gam1).*(nN2out_des.*R.*Temp_E(:,end)).*(((PH./Press_E(:,end)).^gam2)-1) ;
    EnergykWh_des  = (trapz(t5, (Edt_CO2_des + Edt_H2O_des + Edt_N2_des)))/3.6e6 ;
    
% 6) CO2 Compression
    EnergykWh_Comp_vs   = (trapz(t4, (Edt_CO2_vs + Edt_N2_vs)))/3.6e6;
    EnergykWh_Comp_des  = (trapz(t5, (Edt_CO2_des + Edt_N2_des)))/3.6e6 ;
    EnergykWh_Comp     = (EnergykWh_Comp_vs + EnergykWh_Comp_des)/((t4(end)+t5(end))/3600);

% 7) H2O Vacuum Energy 
    EnergykWh_H2O_vs  = (trapz(t4, (Edt_H2O_vs)))/3.6e6 ;
    EnergykWh_H2O_des  = (trapz(t5, (Edt_H2O_des)))/3.6e6 ;
    EnergykW_H2O_vacuum = (EnergykWh_H2O_vs + EnergykWh_H2O_des)/((t4(end)+t5(end))/3600);

 % Total Energy Consumed for a cycle   
    TEnergy_kWh      = (EnergykWh_pz + EnergykWh_cs + EnergykWh_ev + EnergykW_H2O_vacuum + EnergykWh_Comp); %kWh
    TEnergy_tonCO2   = (TEnergy_kWh)/ tonCO2_captured                        ; % Total Energy Consumed per ton (kWh/tonCO2)
    TEnergy_kgCO2   = TEnergy_tonCO2*0.0036  ;
    CO2_Pr          = molCO2_captured/(Vol_Sorbent*(CycleTime/3600))      ;
    H2O_Productivity = molH2O_recovered/(Vol_Sorbent*(CycleTime/3600))      ;
    CO2_Pr_day      = CO2_Pr*((24*0.044)/(1060))  ; % from mol/m3/hr to kgCO2/kgs/day
    CO2kg_s         = (gramCO2_captured/1000)/(CycleTime/3600) ; %kg/hr

% 8) Condenser Work Vapor Stripping
    nH2O_kg = trapz(t4,nH2Oout_vs)*(MM_H2O/1000); %kg/s
   
% Display results in table format
data = {
    'CO2 Sorbed',          sprintf('%.2f g', gramCO2_sorbed);
    'CO2 Recovered',       sprintf('%.2f %%\n'  , CO2Recovery);
    'CO2 Purity',          sprintf('%.2f %%\n'  , CO2purity);
    'H2O Processed',       sprintf('%.2f g', gH2O_used);
    'H2O Lost',            sprintf('%.2f g', gH2O_lost);
    'H2O Recovered',       sprintf('%.2f %%\n' , H2ORecovery);
    'gCO2/gSorbent',       sprintf('%.2f g/g', gCO2_gSorbent);
    'gH2Oused/gCO2',       sprintf('%.2f g/g', gH2Oused_gCO2);
    'gH2Olost/gCO2',       sprintf('%.2f g/g', gH2Olost_gCO2);
    'Electrical Energy',   sprintf('%.2f MJ/kgCO2', TEnergy_kgCO2);
    'CO2 Productivity',    sprintf('%.2f kgCO2/kgs.day', CO2_Pr_day);
    };

% Convert to table and display
PerformanceMetrics = cell2table(data, 'VariableNames', {'Parameter', 'Value'});
fprintf('Output performance metrics:\n');
fprintf('----------------------------------\n');
disp(PerformanceMetrics);
end

%% PLOTS
% % -------------------- Data --------------------
% stages = { ...
%     'CO_2 Sorption\nStep 2', ...
%     'Air Evacuation', ...
%     'Vapor Stripping', ...
%     'Final Desorption', ...
%     'CO_2 Compression'};
% 
% % Build stacked data matrix (Nx2):
% % - Column 1: base component (grey)
% % - Column 2: second component (grey + hatched/shaded on last bar only)
% Y = [ ...
%     ET_ads_MJkgCO2,              0; ...
%     ET_ev_MJkgCO2,               0; ...
%     EnergyMJkgCO2_H2O_vs,        0; ...
%     EnergyMJkgCO2_H2O_des,       0; ...
%     EnergyMJ_comp_vs,  EnergyMJ_comp_des];
% 
% % -------------------- Plot --------------------
% figure(1);   % no clf / clear
% 
% bar_color = [0.6 0.6 0.6];  % Grey for both stacks
% 
% bar = bar(Y, 'stacked', 'FaceColor', 'flat');
% for k = 1:numel(bar)
%     bar(k).CData = repmat(bar_color, size(Y,1), 1);
% end
% 
% % Axes/labels
% ax = gca;
% ax.XTickLabel = stages;
% ax.YAxis.Exponent = 0;
% 
% ylabel({'Electrical Energy'; 'MJ/kgCO_2'}, 'FontSize', 12, 'FontWeight', 'bold');
% xlabel('Process Steps', 'FontSize', 12, 'FontWeight', 'bold');
% title('Electrical Energy Split', 'FontSize', 12, 'FontWeight', 'bold');
% ylim([0 1])
% 
% set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 1.5);
% xtickangle(45);
% 
% % -------------------- Add shaded (hatched) region to one stack --------------------
% % Apply hatch only to the 2nd component on the last bar (CO2 Compression)
% % We draw lines over that segment to "shade" it.
% 
% i = 5; % index of CO2 Compression bar
% x = bar(2).XEndPoints(i);
% w = bar(2).BarWidth;
% 
% y0 = Y(i,1);
% y1 = Y(i,1) + Y(i,2);
% 
% % Hatch settings (adjust spacing for denser/lighter shading)
% dx = w/12;
% x_left  = x - w/2;
% x_right = x + w/2;
% 
% hold on;
% for xx = x_left:dx:x_right
%     plot([xx xx], [y0 y1], 'k-', 'LineWidth', 0.6);  % vertical hatch lines
% end
% hold off;
% 
% % -------------------- Add values on top of each total bar --------------------
% totals = sum(Y, 2);
% for j = 1:numel(totals)
%     text(j, totals(j), sprintf('%.2f', totals(j)), ...
%         'HorizontalAlignment', 'center', ...
%         'VerticalAlignment', 'bottom', ...
%         'FontSize', 10, ...
%         'FontWeight', 'bold');
% end
% 
% % Optional: add a legend for compression split
% legend({'CO_2 comp 1', 'CO_2 comp 2 (shaded)'}, 'Location', 'northwest');


% -------------------- Data --------------------
% stages = { ...
%     'CO_2 Sorption', ...
%     'Air Evacuation', ...
%     'Vapor Stripping', ...
%     'Final Desorption', ...
%     'CO_2 Compression'};
% 
% energy_MJkgCO2 = [ ...
%     ET_ads_MJkgCO2, ...
%     ET_ev_MJkgCO2, ...
%     EnergyMJkgCO2_H2O_vs, ...
%     EnergyMJkgCO2_H2O_des, ...
%     EnergyMJ_CO2comp];
% 
% % Save the energy values as EnCase1
% % save('EnCase1.mat', 'energy_MJkgCO2');
% 
% % -------------------- Plot --------------------
% figure(1);   % no clf, no clearing
% 
% bar_color = [0.6 0.6 0.6];   % Grey
% 
% bb = bar(energy_MJkgCO2, 'FaceColor', 'flat');
% bb.CData = repmat(bar_color, numel(energy_MJkgCO2), 1);
% 
% % Axes settings
% ax = gca;
% ax.XTickLabel = stages;
% ax.YAxis.Exponent = 0;
% 
% ylabel({'Electrical Energy'; 'MJ/kgCO_2'}, 'FontSize', 12, 'FontWeight', 'bold');
% xlabel('Process Steps', 'FontSize', 12, 'FontWeight', 'bold');
% ylim([0 1])
% 
% set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 1.5);
% % xtickangle(45);
% 
% % -------------------- Add values on bars --------------------
% for i = 1:numel(energy_MJkgCO2)
%     text(i, energy_MJkgCO2(i), sprintf('%.2f', energy_MJkgCO2(i)), ...
%         'HorizontalAlignment', 'center', ...
%         'VerticalAlignment', 'bottom', ...
%         'FontSize', 10, ...
%         'FontWeight', 'bold');
% end
