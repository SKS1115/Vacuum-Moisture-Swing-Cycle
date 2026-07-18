
%% Isotherm function to calculate both the Water(H2O) and Carbon Dioxide (CO2) 
% loading on the sorbent, IRA900

function qstar = IsothermMS(YA,YB,P,T,CondParam,isotherm_parameters)

% Condition Parameters (Pre-defined in "ProcessInputParameters" function)
    rho_s  = CondParam(1)  ;              
    R      = CondParam(14) ; 

% Isotherm Parameters (Pre-defined in "ProcessInputParameters" function)
    IEC         = isotherm_parameters(1);
    qm_H2O      = isotherm_parameters(5);
    E_H2O       = isotherm_parameters(6);
    F_H2O       = isotherm_parameters(7);
    G_H2O       = isotherm_parameters(8);
    H_H2O       = isotherm_parameters(9);
   
% The H2O isotherm is based on the temperature Dependent GAB model. 
%  
%                  qm*k*c*x                   E1-E10            E29-E10
%   qH2O =   -----------------------,  c= exp(-------), k= exp(-------)
%            1 - k*x)*(1 + (c-1)*k*x)           RT                RT
%   
%      Emo =-38.6*T + 59338 ; Eml =-50.4*T + 58322 ; Ec = -44.38*T + 58322
%
% where:
% qH2O  : is the concentration of water [mmol/g]
% qm    : is the mono-layer water concentration [mmol/g]
% k_ap  : is an affinity parameter [-]
% c_ap  : is an affinity parameter [-]
% Emo   : is the mono-layer heats of sorption [kJ/mol].
% Eml   : is the heats of sorption of the multi-layer [J/mol]
% Ec    : is the heats of sorption of condensation [kJ/mol]
% T     : is the temperature [K]

% Calculate the water isotherm 
    PH2O     = (YB.*P)                  ;       % Partial Pressure of H2O(Pa)
    PH2O_sat = PressureVaporSat(T)      ;       % Saturation Pressure of H2O at the temperature T
    RH       = PH2O./PH2O_sat           ;       % Relative humidity  
    RH       = min(RH,1)                ;       % Ensure RH doesn't exceed 100%
    Emo      = (-E_H2O.*T) + F_H2O      ;       % Heats of adsorption for 1st layer [J/mol]
    Eml      = (-G_H2O.*T) + H_H2O      ;       % Heats of adsorption for 2nd-9th layer [J/mol]
    Ec       = -44.38.*T + 57220        ;       % Heats of sorption of condensation [J/mol]
    c_ap     = exp((Emo-Ec)./(R.*T))    ;       % affinity parameter [-]
    k_ap     = exp((Eml-Ec)./(R.*T))    ;       % affinity parameter [-]

    qstar2_b  = qm_H2O.*k_ap.*c_ap.*RH                       ;
	qstar2_d  = (1 - (k_ap.*RH)).*(1 + ((c_ap-1).*k_ap.*RH)) ;
    qstar2    = (qstar2_b./qstar2_d)                         ;   % H20 loading on sorbent [mol/kg]  



% The CO2 isotherm is based on a modified langimuir model that has a water dependence with on an exponential factor n. 
% 
%                     IEC*Kms*(H2O^-n)*PCO2
%          qCO2 =   --------------------------
%                    (1 + (Kms*(H2O^-n)*PCO2)        
%   
% where:
% qCO2  : is the CO2 loading at a given water loading [mol/kg]
% IEC   : is the ion exchange capacity or CO2 saturation [mol/kg]
% Kms   : is the equilibrium constant [kg/mol/Pa]
% H2O   : is the concentration of water [mol/kg]
% PCO2  : is the partial pressure of CO2 [Pa]
% n     : is an exponential function fitted from experimental data [-]
    
% Calculate CO2 Isotherm 
% Stephano's fit including 0 point qCO2 (using dynamic temp dependent kms,n)    
    PCO2        = (YA.*P)                       ;  % CO2 partial pressure (Pa)
    kms         = 1e37.*exp(-0.258.*T)          ;
    n           = 1400*exp(-0.02*T)             ;   
    qstar1_b    = (qstar2.^-n).*kms.*PCO2       ;   
    qstar1_d    = 1+((qstar2.^-n).*kms.*PCO2)   ;
    qstar1      = (IEC*qstar1_b)./qstar1_d      ;   % CO2 loading on sorbent mol/kg
    
% Convert CO2 and H2O loading to mol/m3 and export    
    qstar1   = qstar1.*rho_s         ;      % convert CO2 loading on sorbent to [mol/m3];
    qstar2   = qstar2.*rho_s         ;      % convert H2O loading on sorbent to [mol/m3];
    CA       = (YA.*P)./(R.*T)       ;      % CO2 gas concentration [mol/m3] 
    CB       = (YB.*P)./(R.*T)       ;      % H20 gas concentration [mol/m3] 
    
    qstar1(isnan(qstar1)) = 0;
    qstar       = [qstar1, qstar2, RH, CA, CB] ;

end