

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                Energy Balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% load ('ODE_GEN.mat','a','b','c','d','t1','t2','t3','t4','L','N')

%%% FUNCTIONS USED
% j*trapz(i)    : trapz is a numerical integration function that integrates
                % the data set i based on j spacing. 

% Input parameters
    % InputParameters = ProcessInputParameters(N);        % Call input parameter function into this script of code
    % ColumnParam     = InputParameters{1};               % Retrieve respective set of saved parameters from input parameter function
    % CondParam       = InputParameters{2};               % Retrieve respective set of saved parameters from input parameter function
    % isotherm_parameters = InputParameters{3};           % Retrieve respective set of saved parameters from input parameter function

% Retrieve Necessary parameters and calculate constants
    r_in    = ColumnParam(3) ;
    eb      = ColumnParam(5) ;          
    rp      = ColumnParam(7) ;
    rho_s   = CondParam(1)   ;
    C_pg    = CondParam(3)   ;
    C_pa_CO2    = CondParam(4)   ;
    C_ps    = CondParam(5)   ;
    mu      = CondParam(7)   ;
    Dm      = CondParam(8)   ;
    Kz      = CondParam(10)  ; 
    R       = CondParam(14)  ;           
    lambda  = CondParam(16)  ;
    dp      = CondParam(17)  ;
   C_pa_H2O  = CondParam(19);
    PH      = CondParam(20)  ;
    PL      = CondParam(21)  ;
    Pi      = CondParam(22)  ;
    v0      = CondParam(23)  ;
    Ta      = CondParam(25)  ;
    yA0     = CondParam(26)  ;
    yB0     = CondParam(27)  ;
    dUb_CO2 = isotherm_parameters(4) ;
    dUb_H2O  =  isotherm_parameters(10) ;
    
    DL      = 0.7*Dm + ((0.5*v0*dp)/eb) ;
    Viscous = (4/150)*(eb/(1-eb))^2*(rp^2/mu);
    Area    = (pi.*r_in^2)          ;               
    dz      = L/N                   ;
    dz2     = dz/2                  ;
    z       = [0 dz2:dz:L-dz2 L]    ;      % axis from inlet-outlet boundary
    z2      = dz2:dz:L-dz2          ;      % axis only for internal nodes

%% Pressurization
%  Moles IN using ideal gas law 
    Pin_pz     = a(:,2*N+5)                                  ; % Pressure at inlet [Pa]
    Yin_pz     = a(:,1)                                      ; % Mole fraction at inlet is equal to the feed [-]
    Tin_pz     = a(:,5*N+11)                                 ; % Temperature at inlet [K]
    Uin_pz     = -2*Viscous.*((a(:,2*N+6)-a(:,2*N+5))./dz)   ;                                          ; % Velocity at inlet is fixed [m/s] 
    ndot_in_pz = ((Uin_pz.*Pin_pz.*Yin_pz)./Tin_pz)      ; % 
    nCO2_in_pz = (eb*Area./R).*trapz(t1, ndot_in_pz)        ; % [mol]                           

% \Moles OUT using ideal gas law 
    Pout_pz     = a(:,3*N+6)                                  ; % Pressure at outlet [Pa]
    Yout_pz     = a(:,N+2)                                    ; % Mole fraction at outlet [-]
    Tout_pz     = a(:,6*N+12)                                 ; % Temperature at outlet [K]
    Uout_pz     = -2*Viscous.*((a(:,3*N+6)-a(:,3*N+5))./dz)   ; % Velocity at outlet using darcy's law [m/s] 
    Uout_pz     = max(Uout_pz, 0)                            ; % Make negative velocity values = 0. 
    ndot_out_pz = ((Uout_pz.*Pout_pz.*Yout_pz)./Tout_pz)  ; % 
    nCO2_out_pz = (eb*Area./R).*trapz(t1, ndot_out_pz)       ;% [mol]           


%%% Calculate the Heat entering the column 
    rhog_in_pz      = (Pin_pz./Tin_pz).*R^-1                      ; % gas-phase density [mol/m3]
    VTrho_in_pz     = (Uin_pz.*Tin_pz.*rhog_in_pz)              ; % multiply vel*temp*density @inlet
    Heat_in_pz      = eb*Area*C_pg*(trapz(t1, VTrho_in_pz))   % [J/s]

%%% Calculate the Heat leaving the column 
    rhog_out_pz     = (Pout_pz./Tout_pz).*R^-1                    ; % gas-phase density [mol/m3]
    VTrho_out_pz    = (Uout_pz.*Tout_pz.*rhog_out_pz)           ; % multiply vel*temp*density @exit
    Heat_out_pz     = eb*Area*C_pg*(trapz(t1, VTrho_out_pz))  % [J/s]

%%% Calculate the Heat generated        
    dH_qCO2          = -dUb_CO2.*(a(end,3*N+7:4*N+8) - a(1,3*N+7:4*N+8))     ; % dH*qCO2 @t=final - dH*qCO2 @t=initial 
    dH_qH2O           = -dUb_H2O.*(a(end,4*N+9:5*N+10) - a(1,4*N+9:5*N+10))    ; % dH*qH2O @t=final - dH*qH2O @t=initial     
    Heat_gen_pz    = (1-eb)*Area*(trapz(z, dH_qCO2) + trapz(z,dH_qH2O))      

%%% Calculate the Heat in the solid phase 
    Temp_diff_pz   = a(end, 5*N+11:6*N+12) - a(1, 5*N+11:6*N+12)          ; % temp @t=final - temp @t=initial       
    Heat_solid_pz   = (1-eb)*Area*C_ps*rho_s*(trapz(z, Temp_diff_pz))   

%%% Calculate the Heat in the gas of fluid phase
    T_rhog_initial_pz = a(1, 5*N+11:6*N+12).*(a(1,2*N+5:3*N+6)./(R.*a(1, 5*N+11:6*N+12)))       ; % temp*rhog @t=initial 
    T_rhog_final_pz   = a(end, 5*N+11:6*N+12).*(a(end,2*N+5:3*N+6)./(R.*a(end, 5*N+11:6*N+12))) ; % temp*rhog @t=final
    Heat_gas_pz       = eb*Area*C_pg*(trapz(z, (T_rhog_final_pz)-(T_rhog_initial_pz))) 

%%% Calculate the Heat Adsorbed 
    Tq_CO2_initial_pz = a(1,5*N+11:6*N+12).*a(1,3*N+7:4*N+8)       ;  % temp*qCO2 @t=initial 
    Tq_CO2_final_pz   = a(end,5*N+11:6*N+12).*a(end,3*N+7:4*N+8)   ;  % temp*qCO2 @t=final 
    Tq_H2O_initial_pz  = a(1,5*N+11:6*N+12).*a(1,4*N+9:5*N+10)      ;  % temp*qH2O @t=initial 
    Tq_H2O_final_pz    = a(end,5*N+11:6*N+12).*a(end,4*N+9:5*N+10)  ;  % temp*qH2O @t=final 

   Heat_adsorbed_pz  = (1-eb)*Area*(trapz(z, C_pa_CO2.*(Tq_CO2_final_pz-Tq_CO2_initial_pz))...
                       + trapz(z, C_pa_H2O.*(Tq_H2O_final_pz-Tq_H2O_initial_pz)))   

 %%% Calculate the Energy error balance    
    Heat_gain_pz   = Heat_in_pz + Heat_gen_pz                                       ; % total heat energy in and generated
    Heat_lost_pz   = Heat_out_pz + Heat_solid_pz + Heat_gas_pz + Heat_adsorbed_pz ; % total heat energy used and lost
    JENB_pz        = ((abs(Heat_gain_pz - Heat_lost_pz))/Heat_lost_pz)*100          % energy balance error
    
    %% ADSORPTION %%%

% The following lines of code are used to calculate the energy balance for
% the adsorption process. The main notations used are described as follows:
%
% Heat_in_ads       : Heat entering the column [J]
% Heat_out_ads      : Heat leaving the column [J]
% Heat_gen_ads      : Heat generated in the column [J]
% Heat_solid_ads    : Heat accumulated in the solid phase [J]
% Heat_gas_ads      : Heat accumulated in the gas phase [J]
% Heat_adsorbed_ads : Heat adsorbed in the column [J]
% Temp_initial_ads  : Initial column temperature integrated over space [K]
% Temp_final_ads    : Final column temperature integrated over space [K]
% T_rhog_initial_ads  : Initial column gas-phase(GP) density*Temp integrated over space [mol/m3].
% T_rhog_final_ads    : Final column gas-phase(GP) density*Temp integrated over space [mol/m3]

% Major Equations 
% Heat_adsorbed_ads = (1-e)A*C_pa*INT(sum(ncomp)[(T*q_i)final - (T*q_i)initial)dz]
% rhog = pressure/(temp*R)
% T_rhog = temp*rhog = pressure/R


% Adsorption: Moles IN using ideal gas law 
    Pin_ads     = b(:,2*N+5)                                  ; % Pressure at inlet [Pa]
    Yin_ads     = b(:,1)                                      ; % Mole fraction at inlet is equal to the feed [-]
    Tin_ads     = b(:,5*N+11)                                 ; % Temperature at inlet [K]
    Uin_ads     = v0                                          ; % Velocity at inlet is fixed [m/s] 
    ndot_in_ads = ((Uin_ads.*Pin_ads.*Yin_ads)./Tin_ads)      ; % 
    nCO2_in_ads = (eb*Area./R).*trapz(t2, ndot_in_ads)        ; % [mol]                           

% Adsorption: Moles OUT using ideal gas law 
    Pout_ads     = b(:,3*N+6)                                  ; % Pressure at outlet [Pa]
    Yout_ads     = b(:,N+2)                                    ; % Mole fraction at outlet [-]
    Tout_ads     = b(:,6*N+12)                                 ; % Temperature at outlet [K]
    Uout_ads     = -2*Viscous.*((b(:,3*N+6)-b(:,3*N+5))./dz)   ; % Velocity at outlet using darcy's law [m/s] 
    Uout_ads     = max(Uout_ads, 0)                            ; % Make negative velocity values = 0. 
    ndot_out_ads = ((Uout_ads.*Pout_ads.*Yout_ads)./Tout_ads)  ; % 
    nCO2_out_ads = (eb*Area./R).*trapz(t2, ndot_out_ads)       ;% [mol]           


%%% Calculate the Heat entering the column 
    rhog_in_ads      = (Pin_ads./Tin_ads).*R^-1                      ; % gas-phase density [mol/m3]
    VTrho_in_ads     = (Uin_ads.*Tin_ads.*rhog_in_ads)              ; % multiply vel*temp*density @inlet
    Heat_in_ads      = eb*Area*C_pg*(trapz(t2, VTrho_in_ads))   % [J/s]

%%% Calculate the Heat leaving the column 
    rhog_out_ads     = (Pout_ads./Tout_ads).*R^-1                    ; % gas-phase density [mol/m3]
    VTrho_out_ads    = (Uout_ads.*Tout_ads.*rhog_out_ads)           ; % multiply vel*temp*density @exit
    Heat_out_ads     = eb*Area*C_pg*(trapz(t2, VTrho_out_ads))  % [J/s]

%%% Calculate the Heat generated        
    dH_qCO2          = -dUb_CO2.*(b(end,3*N+7:4*N+8) - b(1,3*N+7:4*N+8))     ; % dH*qCO2 @t=final - dH*qCO2 @t=initial 
    dH_qH2O           = -dUb_H2O.*(b(end,4*N+9:5*N+10) - b(1,4*N+9:5*N+10))    ; % dH*qH2O @t=final - dH*qH2O @t=initial     
    Heat_gen_ads     = (1-eb)*Area*(trapz(z, dH_qCO2) + trapz(z,dH_qH2O))      

%%% Calculate the Heat in the solid phase 
    Temp_diff_ads    = b(end, 5*N+11:6*N+12) - b(1, 5*N+11:6*N+12)          ; % temp @t=final - temp @t=initial       
    Heat_solid_ads   = (1-eb)*Area*C_ps*rho_s*(trapz(z, Temp_diff_ads))   

%%% Calculate the Heat in the gas of fluid phase
    T_rhog_initial_ads = b(1, 5*N+11:6*N+12).*(b(1,2*N+5:3*N+6)./(R.*b(1, 5*N+11:6*N+12)))       ; % temp*rhog @t=initial 
    T_rhog_final_ads   = b(end, 5*N+11:6*N+12).*(b(end,2*N+5:3*N+6)./(R.*b(end, 5*N+11:6*N+12))) ; % temp*rhog @t=final
    Heat_gas_ads       = eb*Area*C_pg*(trapz(z, (T_rhog_final_ads)-(T_rhog_initial_ads))) 

%%% Calculate the Heat Adsorbed 
    Tq_CO2_initial_ads = b(1,5*N+11:6*N+12).*b(1,3*N+7:4*N+8)       ;  % temp*qCO2 @t=initial 
    Tq_CO2_final_ads   = b(end,5*N+11:6*N+12).*b(end,3*N+7:4*N+8)   ;  % temp*qCO2 @t=final 
    Tq_H2O_initial_ads  = b(1,5*N+11:6*N+12).*b(1,4*N+9:5*N+10)      ;  % temp*qH2O @t=initial 
    Tq_H2O_final_ads    = b(end,5*N+11:6*N+12).*b(end,4*N+9:5*N+10)  ;  % temp*qH2O @t=final 

    Heat_adsorbed_ads  = (1-eb)*Area*(trapz(z, C_pa_CO2.*(Tq_CO2_final_ads-Tq_CO2_initial_ads))...
                       + trapz(z, C_pa_H2O.*(Tq_H2O_final_ads-Tq_H2O_initial_ads)))   

 %%% Calculate the Energy error balance    
    Heat_gain_ads   = Heat_in_ads + Heat_gen_ads                                       ; % total heat energy in and generated
    Heat_lost_ads   = Heat_out_ads + Heat_solid_ads + Heat_gas_ads + Heat_adsorbed_ads ; % total heat energy used and lost
    JENB_ads        = ((abs(Heat_gain_ads - Heat_lost_ads))/Heat_lost_ads)*100          % energy balance error



%% AIR EVACUTATION %%%

% The following lines of code are used to calculate the energy balance for
% the adsorption process. The main notations used are described as follows:
%
% Heat_in_ads       : Heat entering the column [J]
% Heat_out_ads      : Heat leaving the column [J]
% Heat_gen_ads      : Heat generated in the column [J]
% Heat_solid_ads    : Heat accumulated in the solid phase [J]
% Heat_gas_ads      : Heat accumulated in the gas phase [J]
% Heat_adsorbed_ads : Heat adsorbed in the column [J]
% Temp_initial_ads  : Initial column temperature integrated over space [K]
% Temp_final_ads    : Final column temperature integrated over space [K]
% T_rhog_initial_ads  : Initial column gas-phase(GP) density*Temp integrated over space [mol/m3].
% T_rhog_final_ads    : Final column gas-phase(GP) density*Temp integrated over space [mol/m3]

% Major Equations 
% Heat_adsorbed_ads = (1-e)A*C_pa*INT(sum(ncomp)[(T*q_i)final - (T*q_i)initial)dz]
% rhog = pressure/(temp*R)
% T_rhog = temp*rhog = pressure/R


% Adsorption: Moles IN using ideal gas law 
    Pin_ev     = c(:,2*N+5)                                  ; % Pressure at inlet [Pa]
    Yin_ev     = c(:,N+3)                                      ; % Mole fraction at inlet is equal to the feed [-]
    Tin_ev     = c(:,5*N+11)                                 ; % Temperature at inlet [K]
    Uin_ev     = -2*Viscous.*((c(:,2*N+6)-c(:,2*N+5))./dz)   ; % Velocity at inlet [m/s]                                        ; % Velocity at inlet is fixed [m/s]                          

%  Moles OUT using ideal gas law 
    Pout_ev     = c(:,3*N+6)                                  ; % Pressure at outlet [Pa]
    Yout_ev     = c(:,2*N+4)                                    ; % Mole fraction at outlet [-]
    Tout_ev     = c(:,6*N+12)                                 ; % Temperature at outlet [K]
    Uout_ev     = -2*Viscous.*((c(:,3*N+6)-c(:,3*N+5))./dz)   ; % Velocity at outlet using darcy's law [m/s] 
    Uout_ev     = max(Uout_ev, 0)                            ; % Make negative velocity values = 0. 

%%% Calculate the Heat entering the column 
    rhog_in_ev      = (Pin_ev./Tin_ev).*R^-1                      ; % gas-phase density [mol/m3]
    VTrho_in_ev     = (Uin_ev.*Tin_ev.*rhog_in_ev)              ; % multiply vel*temp*density @inlet
    Heat_in_ev      = eb*Area*C_pg*(trapz(t3, VTrho_in_ev))   % [J/s]

%%% Calculate the Heat leaving the column 
    rhog_out_ev     = (Pout_ev./Tout_ev).*R^-1                    ; % gas-phase density [mol/m3]
    VTrho_out_ev    = (Uout_ev.*Tout_ev.*rhog_out_ev)           ; % multiply vel*temp*density @exit
    Heat_out_ev     = eb*Area*C_pg*(trapz(t3, VTrho_out_ev))  % [J/s]

%%% Calculate the Heat generated        
    dH_qCO2          = -dUb_CO2.*(c(end,3*N+7:4*N+8) - c(1,3*N+7:4*N+8))     ; % dH*qCO2 @t=final - dH*qCO2 @t=initial 
    dH_qH2O           = -dUb_H2O.*(c(end,4*N+9:5*N+10) - c(1,4*N+9:5*N+10))    ; % dH*qH2O @t=final - dH*qH2O @t=initial     
    Heat_gen_ev     = (1-eb)*Area*(trapz(z, dH_qCO2) + trapz(z,dH_qH2O))      

%%% Calculate the Heat in the solid phase 
    Temp_diff_ev    = c(end, 5*N+11:6*N+12) - c(1, 5*N+11:6*N+12)          ; % temp @t=final - temp @t=initial       
    Heat_solid_ev   = (1-eb)*Area*C_ps*rho_s*(trapz(z, Temp_diff_ev))   

%%% Calculate the Heat in the gas of fluid phase
    T_rhog_initial_ev = c(1, 5*N+11:6*N+12).*(c(1,2*N+5:3*N+6)./(R.*c(1, 5*N+11:6*N+12)))       ; % temp*rhog @t=initial 
    T_rhog_final_ev   = c(end, 5*N+11:6*N+12).*(c(end,2*N+5:3*N+6)./(R.*c(end, 5*N+11:6*N+12))) ; % temp*rhog @t=final
    Heat_gas_ev       = eb*Area*C_pg*(trapz(z, (T_rhog_final_ev)-(T_rhog_initial_ev))) 

%%% Calculate the Heat Adsorbed 
    Tq_CO2_initial_ev = c(1,5*N+11:6*N+12).*c(1,3*N+7:4*N+8)       ;  % temp*qCO2 @t=initial 
    Tq_CO2_final_ev   = c(end,5*N+11:6*N+12).*c(end,3*N+7:4*N+8)   ;  % temp*qCO2 @t=final 
    Tq_H2O_initial_ev  = c(1,5*N+11:6*N+12).*c(1,4*N+9:5*N+10)      ;  % temp*qH2O @t=initial 
    Tq_H2O_final_ev    = c(end,5*N+11:6*N+12).*c(end,4*N+9:5*N+10)  ;  % temp*qH2O @t=final 

    Heat_adsorbed_ev  = (1-eb)*Area*(trapz(z, C_pa_CO2.*(Tq_CO2_final_ev-Tq_CO2_initial_ev))...
                       + trapz(z, C_pa_H2O.*(Tq_H2O_final_ev-Tq_H2O_initial_ev)))   

 %%% Calculate the Energy error balance    
    Heat_gain_ev   = Heat_in_ev + Heat_gen_ev                                       ; % total heat energy in and generated
    Heat_lost_ev   = Heat_out_ev + Heat_solid_ev + Heat_gas_ev + Heat_adsorbed_ev ; % total heat energy used and lost
    JENB_ev        = ((abs(Heat_gain_ev- Heat_lost_ev))/Heat_lost_ev)*100           % energy balance error


%% Vapor Stripping 

% The following lines of code are used to calculate the energy balance for
% the adsorption process. The main notations used are described as follows:
%
% Heat_in_ads       : Heat entering the column [J]
% Heat_out_ads      : Heat leaving the column [J]
% Heat_gen_ads      : Heat generated in the column [J]
% Heat_solid_ads    : Heat accumulated in the solid phase [J]
% Heat_gas_ads      : Heat accumulated in the gas phase [J]
% Heat_adsorbed_ads : Heat adsorbed in the column [J]
% Temp_initial_ads  : Initial column temperature integrated over space [K]
% Temp_final_ads    : Final column temperature integrated over space [K]
% T_rhog_initial_ads  : Initial column gas-phase(GP) density*Temp integrated over space [mol/m3].
% T_rhog_final_ads    : Final column gas-phase(GP) density*Temp integrated over space [mol/m3]

% Major Equations 
% Heat_adsorbed_ads = (1-e)A*C_pa*INT(sum(ncomp)[(T*q_i)final - (T*q_i)initial)dz]
% rhog = pressure/(temp*R)
% T_rhog = temp*rhog = pressure/R


% Moles IN using ideal gas law 
    Pin_vp     = e(:,2*N+5)                                  ; % Pressure at inlet [Pa]
    Yin_vp     = e(:,1)                                      ; % Mole fraction at inlet is equal to the feed [-]
    Tin_vp     = e(:,5*N+11)                                 ; % Temperature at inlet [K]
    Uin_vp     = -2*Viscous.*((e(:,2*N+6)-e(:,2*N+5))./dz)                                          ; % Velocity at inlet is fixed [m/s] 
    ndot_in_vp = ((Uin_vp.*Pin_vp.*Yin_vp)./Tin_vp)      ; % 
    nCO2_in_vp = (eb*Area./R).*trapz(t5, ndot_in_vp)         % [mol]                           

% Adsorption: Moles OUT using ideal gas law 
    Pout_vp     = e(:,3*N+6)                                  ; % Pressure at outlet [Pa]
    Yout_vp     = e(:,N+2)                                    ; % Mole fraction at outlet [-]
    Tout_vp     = e(:,6*N+12)                                 ; % Temperature at outlet [K]
    Uout_vp    =  0.1;
    ndot_out_vp = ((Uout_vp.*Pout_vp.*Yout_vp)./Tout_vp)  ; % 
    nCO2_out_vp = (eb*Area./R).*trapz(t5, ndot_out_vp)       % [mol]           


%%% Calculate the Heat entering the column 
    rhog_in_vp      = (Pin_vp./Tin_vp).*R^-1                      ; % gas-phase density [mol/m3]
    VTrho_in_vp     = (Uin_vp.*Tin_vp.*rhog_in_vp)              ; % multiply vel*temp*density @inlet
    Heat_in_vp      = eb*Area*C_pg*(trapz(t5, VTrho_in_vp))   % [J/s]

%%% Calculate the Heat leaving the column 
    rhog_out_vp     = (Pout_vp./Tout_vp).*R^-1                    ; % gas-phase density [mol/m3]
    VTrho_out_vp    = (Uout_vp.*Tout_vp.*rhog_out_vp)           ; % multiply vel*temp*density @exit
    Heat_out_vp     = eb*Area*C_pg*(trapz(t5, VTrho_out_vp))  % [J/s]

%%% Calculate the Heat generated        
    dH_qCO2          = -dUb_CO2.*(e(end,3*N+7:4*N+8) - e(1,3*N+7:4*N+8))     ; % dH*qCO2 @t=final - dH*qCO2 @t=initial 
    dH_qH2O           = -dUb_H2O.*(e(end,4*N+9:5*N+10) - e(1,4*N+9:5*N+10))    ; % dH*qH2O @t=final - dH*qH2O @t=initial     
    Heat_gen_vp     = (1-eb)*Area*(trapz(z, dH_qCO2) + trapz(z,dH_qH2O))      

%%% Calculate the Heat in the solid phase 
    Temp_diff_vp    = e(end, 5*N+11:6*N+12) - e(1, 5*N+11:6*N+12)          ; % temp @t=final - temp @t=initial       
    Heat_solid_vp   = (1-eb)*Area*C_ps*rho_s*(trapz(z, Temp_diff_vp))   

%%% Calculate the Heat in the gas phase
    T_rhog_initial_vp = e(1, 5*N+11:6*N+12).*(e(1,2*N+5:3*N+6)./(R.*e(1, 5*N+11:6*N+12)))       ; % temp*rhog @t=initial 
    T_rhog_final_vp   = e(end, 5*N+11:6*N+12).*(e(end,2*N+5:3*N+6)./(R.*e(end, 5*N+11:6*N+12))) ; % temp*rhog @t=final
    Heat_gas_vp       = eb*Area*C_pg*(trapz(z, (T_rhog_final_vp)-(T_rhog_initial_vp))) 

%%% Calculate the Heat Adsorbed 
    Tq_CO2_initial_vp = e(1,5*N+11:6*N+12).*e(1,3*N+7:4*N+8)       ;  % temp*qCO2 @t=initial 
    Tq_CO2_final_vp   = e(end,5*N+11:6*N+12).*e(end,3*N+7:4*N+8)   ;  % temp*qCO2 @t=final 
    Tq_H2O_initial_vp  = e(1,5*N+11:6*N+12).*e(1,4*N+9:5*N+10)      ;  % temp*qH2O @t=initial 
    Tq_H2O_final_vp    = e(end,5*N+11:6*N+12).*e(end,4*N+9:5*N+10)  ;  % temp*qH2O @t=final 

    Heat_adsorbed_vp  = (1-eb)*Area*(trapz(z, C_pa_CO2.*(Tq_CO2_final_vp-Tq_CO2_initial_vp))...
                       + trapz(z, C_pa_H2O.*(Tq_H2O_final_vp-Tq_H2O_initial_vp)))   

 %%% Calculate the Energy error balance    
    Heat_gain_vp   = Heat_in_vp + Heat_gen_vp                                       ; % total heat energy in and generated
    Heat_lost_vp   = Heat_out_vp + Heat_solid_vp + Heat_gas_vp + Heat_adsorbed_vp ; % total heat energy used and lost
    JENB_vp        = ((abs(Heat_gain_vp - Heat_lost_vp))/Heat_lost_vp)*100           % energy balance error


%% DESORPTION %%%

% The following lines of code are used to calculate the energy balance for
% the desorption process. The main notations used are described as follows:
%
% Heat_in_des       : Heat entering the column [J]
% Heat_out_des      : Heat leaving the column [J]
% Heat_gen_des      : Heat generated in the column [J]
% Heat_solid_des    : Heat accumulated in the solid phase [J]
% Heat_gas_des      : Heat accumulated in the gas phase [J]
% Heat_adsorbed_des : Heat adsorbed in the column [J]
% Temp_initial_des  : Initial column temperature integrated over space [K]
% Temp_final_des    : Final column temperature integrated over space [K]
% rhog_initial_des  : Initial column gas-phase(GP) density integrated over space [mol/m3]
% rhog_final_des    : Final column gas-phase(GP) density integrated over space [mol/m3]

%  Moles IN using ideal gas law 
    Pin_des     = f(:,2*N+5)                                  ; % Pressure at inlet [Pa]
    Yin_des     = f(:,N+3)                                      ; % Mole fraction at inlet is equal to the feed [-]
    Tin_des     = f(:,5*N+11)                                 ; % Temperature at inlet [K]
    Uin_des     = -2*Viscous.*((f(:,2*N+6)-f(:,2*N+5))./dz)   ; % Velocity at inlet [m/s]                                        ; % Velocity at inlet is fixed [m/s]                          

%  Moles OUT using ideal gas law 
    Pout_des     = f(:,3*N+6)                                  ; % Pressure at outlet [Pa]
    Yout_des     = f(:,2*N+4)                                    ; % Mole fraction at outlet [-]
    Tout_des     = f(:,6*N+12)                                 ; % Temperature at outlet [K]
    Uout_des     = -2*Viscous.*((f(:,3*N+6)-f(:,3*N+5))./dz)   ; % Velocity at outlet using darcy's law [m/s] 
    Uout_des     = max(Uout_des, 0)                            ; % Make negative velocity values = 0. 

%%% Calculate the Heat entering the column 
    rhog_in_des      = (Pin_des./Tin_des).*R^-1                     ; % gas-phase density [mol/m3]
    VTrho_in_des     = (Uin_des.*Tin_des.*rhog_in_des)              ; % multiply vel*temp*density @inlet
    Heat_in_des      = eb*Area*C_pg*(trapz(t6, VTrho_in_des))   % [J/s]
    

%%% Calculate the Heat leaving the column 
    rhog_out_des     = (Pout_des./Tout_des).*R^-1                   ; % gas-phase density [mol/m3]
    VTrho_out_des    = (Uout_des.*Tout_des.*rhog_out_des)           ; % multiply vel*temp*density @exit
    Heat_out_des     = eb*Area*C_pg*(trapz(t6, VTrho_out_des))  % [J/s]
    

%%% Calculate the Heat generated        
    dH_qCO2_des      = -dUb_CO2.*(f(end,3*N+7:4*N+8) - f(1,3*N+7:4*N+8))     ; % dH*qCO2 @t=final - dH*qCO2 @t=initial 
    dH_qH2O_des       = -dUb_H2O.*(f(end,4*N+9:5*N+10) - f(1,4*N+9:5*N+10))    ; % dH*qH2O @t=final - dH*qH2O @t=initial     
    Heat_gen_des     = (1-eb)*Area*(trapz(z, dH_qCO2_des) + trapz(z, dH_qH2O_des))     
    

%%% Calculate the Heat in the solid phase 
    Temp_diff_des    = f(end, 5*N+11:6*N+12) - f(1, 5*N+11:6*N+12)          ; % temp @t=final - temp @t=initial       
    Heat_solid_des   = (1-eb)*Area*C_ps*rho_s*(trapz(z, Temp_diff_des))  
            
%%% Calculate the Heat in the gas of fluid phase
    T_rhog_initial_des = f(1, 5*N+11:6*N+12).*(f(1,2*N+5:3*N+6)./(R.*f(1, 5*N+11:6*N+12)))       ; % temp*rhog @t=initial 
    T_rhog_final_des   = f(end, 5*N+11:6*N+12).*(f(end,2*N+5:3*N+6)./(R.*f(end, 5*N+11:6*N+12))) ; % temp*rhog @t=final
    Heat_gas_des       = eb*Area*C_pg*(trapz(z, (T_rhog_final_des)-(T_rhog_initial_des)))   
    
%%% Calculate the Heat Adsorbed 
    Tq_CO2_initial_des = f(1,5*N+11:6*N+12).*f(1,3*N+7:4*N+8)       ;  % temp*qCO2 @t=initial 
    Tq_CO2_final_des   = f(end,5*N+11:6*N+12).*f(end,3*N+7:4*N+8)   ;  % temp*qCO2 @t=final 
    Tq_H2O_initial_des  = f(1,5*N+11:6*N+12).*f(1,4*N+9:5*N+10)      ;  % temp*qH2O @t=initial 
    Tq_H2O_final_des    = f(end,5*N+11:6*N+12).*f(end,4*N+9:5*N+10)  ;  % temp*qH2O @t=final 

    Heat_adsorbed_des  = (1-eb)*Area*(trapz(z, C_pa_CO2.*(Tq_CO2_final_des-Tq_CO2_initial_des))...
                       + trapz(z, C_pa_H2O.*(Tq_H2O_final_des-Tq_H2O_initial_des)))   
       
 %%% Calculate the Energy error balance    
    Heat_gain_des   = Heat_in_des + Heat_gen_des                                   ; % total heat energy in and generated
    Heat_lost_des   = Heat_out_des + Heat_solid_des + Heat_gas_des + Heat_adsorbed_des  ; % total heat energy used and lost
    JENB_des        = abs((Heat_gain_des - Heat_lost_des)/Heat_lost_des)*100           % energy balance error

