

% function PerformanceMetrics = MechanicalEnergy(a, b, c, e, f, t1, t2, t3, t5, t6, RH, N, L, ColumnParam, CondParam)

% Input parameters
    InputParameters     = ProcessInputParameters(N,RH,L);    % Call input parameter function into this script of code
    ColumnParam     = InputParameters{1};               % Retrieve respective set of saved parameters from input parameter function
    CondParam       = InputParameters{2};               % Retrieve respective set of saved parameters from input parameter function
    isotherm_parameters = InputParameters{3};           % Retrieve respective set of saved parameters from input parameter function
    Times              = InputParameters{4};
% Retrieve Necessary parameters and calculate constants
    r_in    = ColumnParam(3) ;
    eb      = ColumnParam(5) ;          
    rp      = ColumnParam(7) ;
    rho_s   = CondParam(1)   ;
    C_pg    = CondParam(3)   ;
    C_pa    = CondParam(4)   ;
    C_ps    = CondParam(5)   ;
    mu      = CondParam(7)   ;
    Dm      = CondParam(8)   ;
    gamma   = CondParam(9)   ;
    Kz      = CondParam(10)  ; 
    R       = CondParam(14)  ;   
    eta     = CondParam(15)  ; 
    lambda  = CondParam(16)  ;
    dp      = CondParam(17)  ;
    Pf      = CondParam(19)  ;
    PH      = CondParam(20)  ;
    PL      = CondParam(21)  ;
    Pi      = CondParam(22)  ;
    v0      = CondParam(23)  ;
    Tf      = CondParam(24)  ;
    Ta      = CondParam(25)  ;
    yA0     = CondParam(26)  ;
    yB0     = CondParam(27)  ;
    dUb_CO2     = isotherm_parameters(4);
    dUb_H20     = isotherm_parameters(10);


    Pvap_in  = PressureVaporSat(Tf) ;
    DL      = 0.7*Dm + ((0.5*v0*dp)/eb) ;
    Viscous = (4/150)*(eb/(1-eb))^2*(rp^2/mu);         
    dz      = L/N                   ;
    dz2     = dz/2                  ;
    z       = [0 dz2:dz:L-dz2 L]    ;      % axis from inlet-outlet boundary

% Calculate Sorbent Properties
    CycleTime   = sum(Times)                 ;   % Total time of cycle
    Area        = eb*pi*(r_in^2)             ;   % Void Area (m2)
    Vol_Sorbent = (1-eb)*pi*(r_in)^2*L       ;   % Volume of sorbent (m3)
    M_CO2       = 44                         ;   % Molar mass of CO2 (g/mol)
    M_H2O       = 18                         ;   % Molar mass of H2O (g/mol)
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
    [yCO2_E, yH2O_E, Press_E, Temp_E, qCO2_E, qH2O_E, yN2_E] = extract_params(e,N);
    [yCO2_F, yH2O_F, Press_F, Temp_F, qCO2_F, qH2O_F, yN2_F] = extract_params(f,N);


%%%.................Extract Key Inputs from Mass Balance Function............................%%%
load('CO2_molarflow.mat');
load('H2O_molarflow.mat');
load('N2_molarflow.mat');
load('CO2_moles');
load('H2O_moles');
load('N2_moles');

% Calculate Key Performance metrics
    molCO2_adsorbed    = molCO2in_pz  + (molCO2in_ads-molCO2out_ads)   ;   % moles of CO2 adsorbed
    molCO2_captured    = molCO2out_vp + molCO2out_des                  ;   % moles of CO2 captured
    gramCO2_captured   = molCO2_captured*M_CO2                         ;   % mass of CO2 captured (grams)
    tonCO2_captured    = gramCO2_captured/1e6                         ;   % mass of CO2 captured (tons)
    gCO2_gSorbent      = gramCO2_captured/g_sorbent                   ;   % grams CO2/grams Sorbent
    CO2purity          = ((molCO2_captured)/(molCO2out_vp + molN2out_vp + molCO2out_des + molN2out_des))*100 ;
    CO2Recovery        = (molCO2_captured/molCO2_adsorbed)*100;
    
    molH2O_used        = molH2Oin_vp                     ;
    molH2O_recovered   = molH2Oout_vp + molH2Oout_des    ; 
    molH2O_lost        = molH2Oout_ads - (molH2Oin_ads + molH2Oin_pz)                   ; 
    H2ORecovery        = (molH2O_recovered/molH2O_used)*100    ; 
    
    gH2O_used          = molH2O_used*M_H2O                 ;  % grams of H2O used 
    gH2O_lost          = molH2O_lost*M_H2O                  ; % grams of H2O lost
    gramsH2O_recovered = molH2O_recovered*M_H2O            ;  % grams of H2O recovered 
    tonsH2O_recovered  = gramsH2O_recovered/1e6            ;  % tons of H2O recovered 
    gH2O_gCO2          = gH2O_used/gramCO2_captured   ;  % grams H2O used/grams CO2 captured
    gH2Oused_gH2Orec   = (gramsH2O_recovered/gH2O_used)*100 ; % grams H2O used/grams H2O recovered
    gH2Olost_gCO2rec = abs(gH2O_lost/gramCO2_captured);

%%%................... Total Energy Consumption...................

%Enery Consumption for Each Step.
    
% 1) Pressurization (Inlet)
    mflowin_pz    = (nCO2in_pz + nH2Oin_pz + nN2in_pz);
    Energydt_pz   = (1/eta)*(gam1).*(mflowin_pz.*R.*Temp_A(:,1)).*(((Press_A(:,1)./PH).^gam2)-1) ;
    EnergykWh_pz  = (trapz(t1, Energydt_pz))/3.6e6    ; % Blower energy in kWh

% 2) Blower Energy during Adsorption 
    mflowin_ads    = (nCO2in_ads + nH2Oin_ads + nN2in_ads);
    Energydt_ads   = ((1/eta)*(gam1).*(mflowin_ads.*R.*Temp_B(:,1)).*(((Press_B(:,1)./PH).^gam2)-1)) ; 
    EnergykWh_ads  = (trapz(t2, Energydt_ads))/3.6e6;     % Blower energy in kWh
    
% 3) Vacuum Energy during Air Evacuation
    mflowout_ev    = (nCO2out_ev + nH2Oout_ev + nN2out_ev);
    Energydt_ev    = (1/eta)*(gam1).*(mflowout_ev.*R.*Temp_C(:,end)).*(((PH./Press_C(:,end)).^gam2)-1) ;
    EnergykWh_ev   = (trapz(t3, Energydt_ev))/3.6e6    ; % Vacuum energy in kWh

% 4) Vacuum Energy during Vapor Stripping
    Edt_CO2_vp     = (1/eta)*(gam1).*(nCO2out_vp.*R.*Temp_E(:,end)).*(((PH./Press_E(:,end)).^gam2)-1) ;
    Edt_H2O_vp     = (1/eta)*(gam1).*(nH2Oout_vp.*R.*Temp_E(:,end)).*(((Pvap_in./Press_E(:,end)).^gam2)-1) ;
    Edt_N2_vp      = (1/eta)*(gam1).*(nN2out_vp.*R.*Temp_E(:,end)).*(((PH./Press_E(:,end)).^gam2)-1) ;
    EnergykWh_vp   = (trapz(t5, (Edt_CO2_vp + Edt_H2O_vp + Edt_N2_vp)))/3.6e6;
    EnergykWh_vp_H2O   = (trapz(t5, (Edt_H2O_vp)))/3.6e6;
  
% 5)  Vacuum Energy during final desorption
    Edt_CO2_des    = (1/eta)*(gam1).*(nCO2out_des.*R.*Temp_F(:,end)).*(((PH./Press_F(:,end)).^gam2)-1) ;
    Edt_H2O_des    = (1/eta)*(gam1).*(nH2Oout_des.*R.*Temp_F(:,end)).*(((Pvap_in./Press_F(:,end)).^gam2)-1) ;
    Edt_N2_des     = (1/eta)*(gam1).*(nN2out_des.*R.*Temp_F(:,end)).*(((PH./Press_F(:,end)).^gam2)-1) ;
    EnergykWh_des  = (trapz(t6, (Edt_CO2_des + Edt_H2O_des + Edt_N2_des)))/3.6e6 ;

 % Total Energy Consumed for a cycle   
    TEnergy_kWh      = (EnergykWh_pz + EnergykWh_ads + EnergykWh_ev + EnergykWh_vp + EnergykWh_des); %kWh
    TEnergy_tonCO2   = (TEnergy_kWh)/ tonCO2_captured                        ; % Total Energy Consumed per ton (kWh/tonCO2)
    TEnergy_kgCO2   = TEnergy_tonCO2*0.0036                        ; % Total Energy Consumed per ton (kWh/tonCO2)
    CO2_Productivity = molCO2_captured/(Vol_Sorbent*(CycleTime/3600))      ;
    CO2_Productivity_day = CO2_Productivity*((24*0.044)/(1060))  ; % from mol/m3/hr to kgCO2/kgs/day
    H2O_Productivity = molH2O_recovered/(Vol_Sorbent*(CycleTime/3600))      ;
        
    % Energy in kWh/tonCO2 for each step
    ET_pz_kWh_tonCO2      = EnergykWh_pz/tonCO2_captured  ;            % Energy pressurization per ton CO2 captured
    ET_ads_kWh_tonCO2     = EnergykWh_ads/tonCO2_captured ;            % Energy adsorption per ton CO2 captured
    ET_ev_kWh_tonCO2      = EnergykWh_ev/tonCO2_captured  ;            % Energy evacuation per ton CO2 captured
    ET_vp_kWh_tonCO2      = EnergykWh_vp/tonCO2_captured ;         % Energy vapor stripping per ton CO2 captured
    ET_des_kWh_tonCO2     = EnergykWh_des/tonCO2_captured ;            % Energy desorption per ton CO2 captured
       

  save ("Step_Energy_kWhtn.mat",'ET_pz_kWh_tonCO2','ET_ads_kWh_tonCO2','ET_ev_kWh_tonCO2','ET_vp_kWh_tonCO2','ET_des_kWh_tonCO2','-v7.3')

% Display results in table format
data = {
    'CO2 Adsorbed',          sprintf('%.2f mol', molCO2_adsorbed);
    'CO2 Recovered',         sprintf('%.2f %%\n'  , CO2Recovery);
    'CO2 Purity',            sprintf('%.2f %%\n'  , CO2purity);
    'H2O Used',              sprintf('%.2f g', gH2O_used);
    'H2O Lost',              sprintf('%.2f g', gH2O_lost);
    'H2O Recovered',         sprintf('%.2f %%\n' , H2ORecovery);
    'gCO2cap/gSorbent',         sprintf('%.2f g/g', gCO2_gSorbent);
    'gH2Oused/gCO2cap',             sprintf('%.2f g/g', gH2O_gCO2);
    'gH2Olost/gCO2cap',      sprintf('%.2f g/g', gH2Olost_gCO2rec);
    'Energy Consumed',       sprintf('%.2f MJ/kgCO2', TEnergy_kgCO2);
    'CO2 Productivity',      sprintf('%.2f kgCO2/kgs.day', CO2_Productivity_day);
    };

% Convert to table and display
PerformanceMetrics = cell2table(data, 'VariableNames', {'Parameter', 'Value'});
disp(PerformanceMetrics);


% Plot
% Define the stages and corresponding energy values (in MWhr/ton of CO2)
% load Step_Energy_kWhtn.mat

% Define the stages and corresponding energy values (in MWhr/ton of CO2)
stages = {'Pressurization', 'CO2 Adsorption', 'Air Evacuation', 'Vapor Stripping', 'Final Desorption'};
energyt_values_kWhCO2 = [ET_pz_kWh_tonCO2, ET_ads_kWh_tonCO2, ET_ev_kWh_tonCO2, ET_vp_kWh_tonCO2, ET_des_kWh_tonCO2];

% Define colors for the groups
blower_color = [0.3, 0.6, 0.9];  % Light blue for blower stages
vacuum_color = [0.8, 0.4, 0.4];  % Light red for vacuum stages

%........... Create the bar graph with custom colors for kWh/tonCO2.......
figure(1);
graph = bar(energyt_values_kWhCO2, 'FaceColor', 'flat');
graph.CData(1:2, :) = repmat(blower_color, 2, 1);   % Apply blower color to Pressurization and Adsorption
graph.CData(3:5, :) = repmat(vacuum_color, 3, 1);   % Apply vacuum color to Evacuation, Evaporation, and Desorption

% Set the x-axis labels & Label the axes and title
set(gca, 'XTickLabel', stages); ax = gca;   ax.YAxis.Exponent = 0;
ylabel('kWhr/ton CO_2');    xlabel('Process Stages');
title('Electrical Energy');
% ylim([0 500])

% Add custom legend by plotting invisible bars with appropriate colors
hold on;
bar(nan, 'FaceColor', blower_color);  % Invisible bar for blower legend
bar(nan, 'FaceColor', vacuum_color);  % Invisible bar for vacuum legend
hold off;
legend({'Blower', '','Vacuum'}, 'Location', 'northeast');
xtickangle(45);

% Add values on top of each bar
for i = 1:length(energyt_values_kWhCO2)
    text(i, energyt_values_kWhCO2(i) + 1, sprintf('%.1f', energyt_values_kWhCO2(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
end
hold off;

%........... Create the bar graph with custom colors for MJ/kgCO2.......
figure(2);
energyt_values_MJkgCO2 = energyt_values_kWhCO2*0.0036;
graph = bar(energyt_values_MJkgCO2, 'FaceColor', 'flat');
graph.CData(1:2, :) = repmat(blower_color, 2, 1);   % Apply blower color to Pressurization and Adsorption
graph.CData(3:5, :) = repmat(vacuum_color, 3, 1);   % Apply vacuum color to Evacuation, Evaporation, and Desorption

% Set the x-axis labels & Label the axes and title
set(gca, 'XTickLabel', stages); 
ax = gca;   
ax.YAxis.Exponent = 0;

ylabel('MJ/kgCO_2', 'FontSize', 12, 'FontWeight', 'bold');    
xlabel('Process Stages', 'FontSize', 12, 'FontWeight', 'bold');
title('Electrical Energy', 'FontSize', 12, 'FontWeight', 'bold');
% ylim([0 2])

% Make axis tick labels bold and boundaries thick
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 1.5);

% Add custom legend by plotting invisible bars with appropriate colors
hold on;
bar(nan, 'FaceColor', blower_color);  % Invisible bar for blower legend
bar(nan, 'FaceColor', vacuum_color);  % Invisible bar for vacuum legend
hold off;
legend({'Blower', '','Vacuum'}, 'Location', 'northeast');

% Rotate x-tick labels
xtickangle(45);

% Add values on top of each bar
for i = 1:length(energyt_values_MJkgCO2)
    text(i, energyt_values_MJkgCO2(i), sprintf('%.1f', energyt_values_MJkgCO2(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 10, 'FontWeight', 'bold');
end
hold off;
