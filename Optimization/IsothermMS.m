
% Isotherm function to calculate both the Water(H2O) and Carbon Dioxide (CO2) Isotherms for Moisture Swing Material, IRA900

function qstar = IsothermMS(YA,YB,P,T,CondParam,isotherm_parameters)

% Condition Parameters (Pre-defined in "ProcessInputParameters" function)
    rho_s  = CondParam(1)  ;              
    R      = CondParam(14) ; 

% Isotherm Parameters (Pre-defined in "ProcessInputParameters" function)
    IEC         = isotherm_parameters(1);
    % kms         = isotherm_parameters(2);
    % n           = isotherm_parameters(3);
    qm_H2O      = isotherm_parameters(5);
    C_H2O       = isotherm_parameters(6);
    D_H2O       = isotherm_parameters(7);
    F_H2O       = isotherm_parameters(8);
    G_H2O       = isotherm_parameters(9);
   
% The H2O isotherm is based on the temperature Dependent GAB model. 
%  
%                  qm*k*c*x                   E1-E10            E29-E10
%   qH2O =   -----------------------,  c= exp(-------), k= exp(-------)
%            1 - k*x)*(1 + (c-1)*k*x)           RT                RT
%   
%      E1 =-38.6*T + 59338 ; E29 =-50.4*T + 58322 ; E10 = -44.38*T + 58322
%
% where:
% qH2O  : is the cocentration of water [mmol/g]
% qm    : is the monolayer water concentration [mmol/g]
% k     : is the equilibrium constant [mol/kg*Pa]
% c     : is the concentration of water [mmol/g]
% x     : is the partial pressure of CO2 [Pa]
% n     : is an exponential function fitted from experimental data [-] 
% E1    : is the monolayer heats of adsorption [kJ/mol].
% E29   : is the heats of adsorption of the 2nd to 9th layer [kJ/mol]
% E10   : is the heats of adsorption of the 10th and higher layers [kJ/mol]
% T     : is the temperature [K]

% Calculate the water isotherm 
    PH2O     = (YB.*P)                  ;       % Partial Pressure of H2O(Pa)
    PH2O_sat = PressureVaporSat(T)      ;       % Saturation Pressure of H2O at a given T
    RH       = PH2O./PH2O_sat           ;       % Relative humidity  
    RH       = min(RH,1);
    E1       = (-C_H2O.*T) + D_H2O      ;       % Heats of adsorption for 1st layer [J/mol]
    E29      = (-F_H2O.*T) + G_H2O      ;       % Heats of adsorption for 2nd-9th layer [J/mol]
    E10      = -44.38.*T + 57220        ;       % Heats of adsorption for 10th+ layer [J/mol]
    c_ap     = exp((E1-E10)./(R.*T))    ;       % affinity parameter [-]
    k_ap     = exp((E29-E10)./(R.*T))   ;       % affinity parameter [-]

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
% qCO2  : is the CO2 loading at any given water loading [mmol/g]
% IEC   : is the ion exchange capacity or CO2 saturation [mmol/g]
% Kms   : is the equilibrium constant [mol/kg*Pa]
% H2O   : is the concentration of water [mmol/g]
% PCO2  : is the partial pressure of CO2 [Pa]
% n     : is an exponential function fitted from experimental data [-]
    
% Calculate CO2 Isotherm 
    
    PCO2        = (YA.*P)                                   ;      % CO2 partial pressure (Pa)

 % Stephano's fit including 0 point qCO2
    kms         = 1e37.*exp(-0.258.*T);
    n           = 1400*exp(-0.02*T);   
    

    qstar1_b    = (qstar2.^-n).*kms.*PCO2                   ;   
    qstar1_d    = 1+((qstar2.^-n).*kms.*PCO2)               ;
    qstar1      = (IEC*qstar1_b)./qstar1_d                  ;       
    
    CA       = (YA.*P)./(R.*T)       ;      % CO2 concentration [mol/m3] 
    CB       = (YB.*P)./(R.*T)       ;      % H20 concentration [mol/m3] 

    
% Convert CO2 and H2O concentrations to mol/m3 and export    
    qstar1      = qstar1.*rho_s              ;       % CO2 loading on sorbent [mol/m3];
    qstar2      = qstar2.*rho_s              ;       % H2O loading on sorbent [mol/m3]; 
    
    qstar1(isnan(qstar1)) = 0;
    qstar       = [qstar1, qstar2, RH, CA, CB] ;

end