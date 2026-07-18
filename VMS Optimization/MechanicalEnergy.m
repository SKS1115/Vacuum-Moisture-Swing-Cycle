

function Objectives = MechanicalEnergy(a, b, c, d, e, t1, t2, t3, t4, t5, RH, N, L, ColumnParam, CondParam, Times)

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

% Calculate Key Performance metrics
    molCO2_adsorbed    = molCO2in_pz  + (molCO2in_cs-molCO2out_cs)   ;   % moles of CO2 adsorbed
    molCO2_captured    = molCO2out_vs + molCO2out_des                ;   % moles of CO2 captured
    gramCO2_captured   = molCO2_captured*M_CO2                       ;   % mass of CO2 captured (grams)
    tonCO2_captured    = gramCO2_captured/1e6                        ;   % mass of CO2 captured (tons)
    gCO2_gSorbent      = gramCO2_captured/g_sorbent                  ;   % grams CO2/grams Sorbent
    CO2purity          = ((molCO2_captured)/(molCO2out_vs + molN2out_vs + molCO2out_des + molN2out_des))*100 ;
    CO2Recovery        = (molCO2_captured/molCO2_adsorbed)*100;

    molH2O_used        = molH2Oin_vs                     ;
    molH2O_recovered   = molH2Oout_vs + molH2Oout_des    ; 
    molH2O_lost        = molH2Oout_cs - (molH2Oin_cs + molH2Oin_pz) ;  
    H2ORecovery        = (molH2O_recovered/molH2O_used)*100    ; 
    
    gramsH2O_used      = molH2O_used*M_H2O                 ;  % grams of H2O used 
    gramsH2O_recovered = molH2O_recovered*M_H2O            ;  % grams of H2O recovered 
    gramsH2O_lost      = molH2O_lost*M_H2O                 ;  % grams of H2O lost during adsorption
    tonsH2O_recovered  = gramsH2O_recovered/1e6            ;  % tons of H2O recovered 
    gH2Oused_gCO2cap   = gramsH2O_used/gramCO2_captured    ;  % grams H2O recovered/grams CO2 captured
    gH2Olost_gCO2cap   = gramsH2O_lost/gramCO2_captured    ;
 
% Early termination for nonphysical cases
    if gH2Olost_gCO2cap < 0 || isnan(gH2Olost_gCO2cap) || isinf(gH2Olost_gCO2cap)
        MECH_ENERGY        = 1e6;   % Large energy penalty
        CO2_PRODUCTIVITY   = 0  ;   % No productivity
        gH2Olost_gCO2cap   = 1e6;   % Large water loss penalty
        gH2Oused_gCO2cap   = 1e6;   % Large penalty
        CO2_Purity         = 0  ;   % Set purity to 0

        fprintf('⚠️  Skipping case: Negative/invalid water loss detected.\n');
        Objectives = [MECH_ENERGY , CO2_PRODUCTIVITY, gH2Olost_gCO2cap, gH2Oused_gCO2cap, CO2_Purity] ; 

        return;  % Exit this function immediately
    end

%%%................... Total Energy Consumption...................

%Enery Consumption for Each Step.
    
% 1) Pressurization (Inlet)
    mflowin_pz    = (nCO2in_pz + nH2Oin_pz + nN2in_pz);
    Energydt_pz   = (1/eta)*(gam1).*(mflowin_pz.*R.*Temp_A(:,1)).*(((Press_A(:,1)./PH).^gam2)-1) ;
    EnergykWh_pz  = (trapz(t1, Energydt_pz))/3.6e6    ; % Blower energy in kWh

% 2) Blower Energy during Adsorption 
    mflowin_cs    = (nCO2in_cs + nH2Oin_cs + nN2in_cs);
    Energydt_cs   = ((1/eta)*(gam1).*(mflowin_cs.*R.*Temp_B(:,1)).*(((Press_B(:,1)./PH).^gam2)-1)) ; 
    EnergykWh_cs  = (trapz(t2, Energydt_cs))/3.6e6    ; % Blower energy in kWh
    
% 3) Vacuum Energy during Air Evacuation
    mflowout_ev    = (nCO2out_ev + nH2Oout_ev + nN2out_ev);
    Energydt_ev    = (1/eta)*(gam1).*(mflowout_ev.*R.*Temp_C(:,end)).*(((PH./Press_C(:,end)).^gam2)-1) ;
    EnergykWh_ev   = (trapz(t3, Energydt_ev))/3.6e6    ; % Vacuum energy in kWh

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

 % Total Energy Consumed for a cycle   
    TEnergy_kWh      = (EnergykWh_pz + EnergykWh_cs + EnergykWh_ev + EnergykWh_vs + EnergykWh_des); %kWh
    TEnergy_tonCO2   = (TEnergy_kWh)/ tonCO2_captured                        ; % Total Energy Consumed per ton (kWh/tonCO2)
    CO2_Productivity = molCO2_captured/(Vol_Sorbent*(CycleTime/3600))      ;
    H2O_Productivity = molH2O_recovered/(Vol_Sorbent*(CycleTime/3600))      ;
        
    % Energy in kWh/tonCO2 for each step
    ET_pz_kWh_tonCO2    = EnergykWh_pz/tonCO2_captured  ;         % Energy pressurization per ton CO2 captured
    ET_cs_kWh_tonCO2    = EnergykWh_cs/tonCO2_captured  ;         % Energy adsorption per ton CO2 captured
    ET_ev_kWh_tonCO2    = EnergykWh_ev/tonCO2_captured  ;         % Energy evacuation per ton CO2 captured
    ET_vs_kWh_tonCO2    = EnergykWh_vs/tonCO2_captured  ;         % Energy vapor stripping per ton CO2 captured
    ET_des_kWh_tonCO2   = EnergykWh_des/tonCO2_captured ;         % Energy desorption per ton CO2 captured
       
% Display results in table format
data = {
    'CO2 Adsorbed',          sprintf('%.4f mol', molCO2_adsorbed);
    'CO2 Recovered',         sprintf('%.4f %%\n'  , CO2Recovery);
    'CO2 Purity',            sprintf('%.4f %%\n'  , CO2purity);
    'H2O Used',              sprintf('%.4f mol', molH2O_used);
    'H2O Recovered',         sprintf('%.4f %%\n' , H2ORecovery);
    'gCO2/gSorbent',         sprintf('%.4f g/g', gCO2_gSorbent);
    'gH2O/gCO2',             sprintf('%.4f g/g', gH2Oused_gCO2cap);
    'molH2Olost/molCO2rec',  sprintf('%.4f g/g', gH2Olost_gCO2cap);
    'Energy Consumed',       sprintf('%.4f kWh/tonCO2', TEnergy_tonCO2);
    'CO2 Productivity',      sprintf('%.4f molCO2/m3.hr', CO2_Productivity);
    'H2O Productivity',      sprintf('%.4f molH2O/m3.hr', H2O_Productivity); 
    };

% Convert to table and display
PerformanceMetrics = cell2table(data, 'VariableNames', {'Parameter', 'Value'});

% Define the stages and corresponding energy values (in MWhr/ton of CO2)
stages = {'Pressurization', 'CO2 Adsorption', 'Air Evacuation', 'Vapor Stripping', 'Final Desorption'};
energyt_values_kWhCO2 = [ET_pz_kWh_tonCO2, ET_cs_kWh_tonCO2, ET_ev_kWh_tonCO2, ET_vs_kWh_tonCO2, ET_des_kWh_tonCO2];

% Define colors for the groups
blower_color = [0.3, 0.6, 0.9];  % Light blue for blower stages
vacuum_color = [0.8, 0.4, 0.4];  % Light red for vacuum stages


Objectives = [TEnergy_tonCO2 , CO2_Productivity, gH2Olost_gCO2cap, gH2Oused_gCO2cap, CO2purity] ; 

end