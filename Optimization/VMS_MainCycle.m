
function [MECH_ENERGY, CO2_PRODUCTIVITY, gH2Olost_gCO2cap, gH2Oused_gCO2cap, CO2_Purity] = VMS_MainCycle(N, RH, t_ad, t_vp, t_de, v_ads, v_vp)
tic 

L = 0.01;
dz = L/N;
RH2 = RH;

% Input parameters
    InputParameters = ProcessInputParameters(N,RH,t_ad,t_vp,t_de,v_ads, v_vp); % Call input parameter function into this script of code
    ColumnParam         = InputParameters{1};               % Retrieve respective set of saved parameters from input parameter function
    CondParam           = InputParameters{2};               % Retrieve respective set of saved parameters from input parameter function
    isotherm_parameters = InputParameters{3};           % Retrieve respective set of saved parameters from input parameter function
    Times               = InputParameters{4};               % Retrieve respective set of saved parameters from input parameter function

% Retrieve Necessary parameters and calculate constants
    eb      = ColumnParam(5) ;
    rp      = ColumnParam(7) ;
    C_pg    = CondParam(3)   ;
    mu      = CondParam(7)   ;
    Dm      = CondParam(8)   ;
    Kz      = CondParam(10)  ; 
    R       = CondParam(14)  ;
    lambda  = CondParam(16)  ;
    dp      = CondParam(17)  ;
    deltaP  = CondParam(19)  ;
    PH      = CondParam(20)  ;
    PL      = CondParam(21)  ;
    Pi      = CondParam(22)  ;
    v0      = CondParam(23)  ;
    Tf      = CondParam(24)  ;
    Ta      = CondParam(25)  ;
    yA0     = CondParam(26)  ;
    yB0     = CondParam(27)  ;
    yC0     = CondParam(28)  ;
    u_vp    = CondParam(29)  ;
    MW_air   = 0.02896;
    MW_water = 0.018015; 
    DL      = 0.7*Dm;    
    viscous_term  = (150*(1-eb)^2*mu)/(eb^2*dp^2) ; 
    dz2     = dz/2           ;


% Call PSA Cycle Functions
    Pressurization_fxn    = @(t, x) FuncPressurization(t, x,dz,N,L,ColumnParam,CondParam,isotherm_parameters) ;
    Adsorption_fxn        = @(t, x) FuncAdsorption(t, x, dz,N,L,ColumnParam,CondParam,isotherm_parameters) ;
    Evacuation_fxn        = @(t, x) FuncEvacuation(t, x, dz,N,L,ColumnParam,CondParam,isotherm_parameters) ;
    Desorption_fxn        = @(t, x) FuncDesorption(t, x, dz,N,L,ColumnParam,CondParam,isotherm_parameters) ;

% Retrieve PSA Cycle Times
    t_pz    = Times(1);                       
    t_ads   = Times(2);                      
    t_evac  = Times(3); 
    t_evap1 = Times(4); 
    t_evap2 = Times(5); 
    t_des   = Times(6);  

    
% Initialize the column 
    qstar               = IsothermMS(yA0, yB0, PL, Ta, CondParam, isotherm_parameters)  ; %(yCO2,yH20,Press, Temp, isotherm input)
    x0                  = zeros(7*N+14,1) ;                 % Initialize the size of the IC matrix
    x0(1:N+2)           = yA0               ;                   % Intial conditions for CO2 mole fraction
    x0(N+3:2*N+4)       = yB0             ;                 % Intial conditions for H20 mole fraction
    x0(2*N+5:3*N+6)     = PL              ;                 % Intial conditions for total pressure
    x0(3*N+7:4*N+8)     = qstar(1)        ;                 % Intial conditions for CO2 reaction rate
    x0(4*N+9:5*N+10)    = qstar(2)        ;                 % Intial conditions for H20 reaction rate
    x0(5*N+11:6*N+12)   = Ta              ;                 % Intial conditions for total temperature
    x0(6*N+13:7*N+14)   = Ta              ;   
    
 %% Begin simulating aVMS cycle.
%   Run the simulation until the change in the temperature, gas mole fraction
%   and CO2 molar loading is less than 0.5%. For the molar loading and the
%   mole fraction, the simulation also stops if the absolute change in the
%   state variable is less than 5e-5 and 2.5e-4 respectively. 

    % opts1 = odeset( 'RelTol', 1e-3) ;
    % opts2 = odeset( 'RelTol', 1e-3) ;
    % opts3 = odeset( 'RelTol', 1e-3) ;
    % opts4 = odeset( 'RelTol', 1e-3) ;
    % opts5 = odeset( 'RelTol', 1e-3) ;

    opts1 = odeset('RelTol',1e-6, 'AbsTol',1e-9);
    opts2 = odeset('RelTol',1e-6, 'AbsTol',1e-9); 
    opts3 = odeset('RelTol',1e-6, 'AbsTol',1e-9); 
    opts4 = odeset('RelTol',1e-6, 'AbsTol',1e-9); 
    opts5 = odeset('RelTol',1e-6, 'AbsTol',1e-9); 

    
cycles = 30;
cycle_data   = zeros(cycles,8) ; % This will hold the cycle number, CSS, and mass balance
CSS_tolerance = 0.005;
MBS_tolerance = 0.005;
break_condition_met = false    ; % Track whether the break condition is et

for i = 1:cycles

     cycle_start = tic; % Track time for each cycle

 %%% 1. SIMULATE PRESSURIZATION STEP [OPEN--CLOSED]
        [t1, a] = ode15s(Pressurization_fxn, [0 t_pz], x0, opts1) ;

        % Correct the output to match boundary conditions (clean up results from simulation)
        %%% Inlet 
        idx             = find(a(:, 2*N+5) < a(:,2*N+5))    ; % Find values of P(1) < P(2)
        a(idx, 2*N+5)   = a(idx, 2*N+6)                    ; 
        MW_in           = MW_air + (MW_water-MW_air).*a(:, N+3);
        rho_g1     = a(:,2*N+5)/(R*Ta);  
        kinetic_term  = (1.75*(1-eb).*rho_g1.*MW_in)./(eb*dp) ;
        dsmt1   = viscous_term^2 - (4*kinetic_term.*(a(:,2*N+6)-a(:,2*N+5))*(2/dz)) ;
        u    = (-viscous_term + sqrt(dsmt1))./(2*kinetic_term) ;     
        a(:, 1)         = (a(:,2) + yA0.*((u.*dz)./(2*DL)))./ (1 + ((u.*dz)./(2*DL)))            ; 
        a(:, N+3)       = (a(:,N+4) + yB0.*((u.*dz)./(2*DL)))./ (1 + ((u.*dz)./(2*DL)))            ; %           
        rho_g           = a(:, 2*N+5)./(R*Ta)    ; 
        Constant2       = rho_g.*((eb.*u.*C_pg*dz)./(2*Kz)) ;
        a(:, 5*N+11)    = (a(:, 5*N+12) + (Ta.*Constant2(:)))./(1+Constant2(:));      % T(1) (Dankwerts BC)
        a(:, 3*N+7)     = a(:, 3*N+8)                       ; % qCO2(1) = qCO2(2)
        a(:, 4*N+9)     = a(:, 4*N+10)                      ; % qH20(1)  = H20(2)

        %%% Outlet
        a(:, N+2)       = a(:, N+1)                         ; % outlet yCO2(N+2) = yCO2(N+1), imposing the zero gradient BC
        a(:, 2*N+4)     = a(:, 2*N+3)                       ; % yH20(N+2) = yH20(N+1)
        a(: ,3*N+6)     = a(: ,3*N+5)                       ; % P
        a(:, 4*N+8)     = a(:, 4*N+7)                       ; % qCO2(N+2) = qCO2(N+1)
        a(:, 5*N+10)    = a(:, 5*N+9)                       ; % qH20(N+2)  = qH20(N+1)        
        a(:, 6*N+12)    = a(:, 6*N+11)                      ; % T(N+2) = T(N+1) 
        yCa = 1-(a(:,1:N+2) + a(:,N+3:2*N+4))               ; % yN2(1:N+2)
       
        % Prepare initial conditions for Adsorption Step
        x10         = a(end, :)'  ; % Final state of previous step is the initial state for current step
        % x10(1)      = yA0;
        % x10(N+3)    = yB0;
        x10(2*N+5:3*N+6)  = PH;
        % x10(5*N+11) = Ta;

        % Initial conditions of states at first step of the PSA Cycle
        statesIC = a(1, [2:N+1, N+4:2*N+3, 2*N+6:3*N+5, 3*N+8:4*N+7, 4*N+10:5*N+9, 5*N+12:6*N+11]);  

  %%% 2. SIMULATE ADSORPTION STEP [OPEN--OPEN]
        [t2, b] = ode15s(Adsorption_fxn, [0 t_ads], x10,opts2) ; % When initial conditions are based on pressurization
        
        % Correct the output (clean up results from simulation)
        %%% Inlet        
        b(:, 5*N+11)    = Ta; 
        b(:, 1)         = (b(:,2) + yA0*((v0*dz)/(2*DL)))/ (1 + ((v0*dz)/(2*DL)))            ; % yCO2(1)   = dankwerts BC
        b(:, N+3)       = (b(:,N+4) + yB0*((v0*dz)/(2*DL)))/ (1 + ((v0*dz)/(2*DL)))            ; % yH20(1)   = dankwerts BC
        % rho_g           = PH/(R*Ta)    ; % P(1)/(R*T(1))
        % Constant2       = rho_g*((eb*v0*C_pg*dz)/(2*Kz)) ;
        % b(:, 5*N+11)    = (b(:, 5*N+12) + (Ta.*Constant2(:)))./(1+Constant2(:));      % T(1) (Dankwerts BC     
        
        
        MW_in           = MW_air + (MW_water-MW_air)*b(:, N+3);
        b(:, 2*N+5)    = ((viscous_term*v0*dz2) + b(:, 2*N+6))/(1-(dz/2/R/b(:, 5*N+11))*(1.75*(1-eb)*MW_in/dp/eb)*v0^2);  
        b(:, 3*N+7)     = b(:, 3*N+8)                       ; % qCO2(1)   = qCO2(2)
        b(:, 4*N+9)     = b(:, 4*N+10)                      ; % qH20(1)    = qH20(2)

        %%% Outlet
        b(:, 3*N+6)     = PH;   
        idx             = find(b(:, 3*N+5) < PH)            ; % find P(N+1) <PH = P(N+2)
        b(idx, 3*N+6)   = b(idx, 3*N+5)                     ; % P(N+2) = P(N+1)
        b(:, N+2)       = b(:, N+1)                         ; % yCO2(N+2) = yCO2(N+1)
        b(:, 2*N+4)     = b(:, 2*N+3)                       ; % yH20(N+2) = yCO2(N+1)     
        b(:, 4*N+8)     = b(:, 4*N+7)                       ; % qCO2(N+2) = qCO2(N+1)
        b(:, 5*N+10)    = b(:, 5*N+9)                       ; % qH20(N+2)  = qH20(N+1)        
        b(:, 6*N+12)    = b(:, 6*N+11)                      ; % T(N+2) = T(N+1)                   
        yCb             = 1-(b(:,1:N+2) + b(:,N+3:2*N+4))   ; % yN2(1:N+2)
              
        % Prepare initial conditions for N2 Evacuation Step
        x20         = b(end, :)'  ; % Final state of previous step is the initial state for current step
        
     
  %%% 3. SIMULATE AIR EVACUATION STEP [CLOSED--OPEN]
        [t3, c] = ode15s(Evacuation_fxn, [0 t_evac], x20, opts3) ;
     
        % Correct the Boundary Conditions (clean up results from simulation)
        %%% Inlet
        c(:, 2*N+5)     = c(:, 2*N+6)                       ; % P(1)   = P(2)
        c(:, 1)         = c(:, 2)                           ; % yCO2(1) = yCO2(2)
        c(:, N+3)       = c(:, N+4)                         ; % yH20(1) = yH20(2)
        c(:, 3*N+7)     = c(:, 3*N+8)                       ; % qCO2(1)   = qCO2(2)
        c(:, 4*N+9)     = c(:, 4*N+10)                      ; % qH20(1)    = qH20(2)
        c(:, 5*N+11)    = c(:, 5*N+12)                      ; % T(1) = T(2)
       
        %%% Outlet
        idx             = find(c(:, 3*N+6) > c(:,3*N+5))    ; % Find values of P(N+1) < P(N+2)
        c(idx, 3*N+6)   = c(idx, 3*N+5)                     ; % P(N+2)    = P(N+1) for P(N+1) values < P(N+2)
        c(:, N+2)       = c(:, N+1)                         ; % yCO2(N+2) = yCO2(N+1)
        c(:, 2*N+4)     = c(:, 2*N+3)                       ; % yH20(N+2) = yH20(N+1)       
        c(:, 4*N+8)     = c(:, 4*N+7)                       ; % qCO2(N+2) = qCO2(N+1)
        c(:, 5*N+10)    = c(:, 5*N+9)                       ; % qH20(N+2)  = qH20(N+1)
        c(:, 6*N+12)    = c(:, 6*N+11)                      ; % T(N+2) = T(N+1)
        yCc = 1-(c(:,1:N+2) + c(:,N+3:2*N+4))               ; % yN2(1:N+2)

        % Prepare initial conditions for Closed-Closed
        x30             = c(end, :)'                ; % Final state of previous step is the initial state for current step        

        % Set mole fraction of water based on RH
        RH              = 100                  ;
        Plast           = x30(3*N+6) ; 
        Pvap_in         = PressureVaporSat(Tf) ;
        yB0             = (RH*Pvap_in)/(100*Pvap_in)      ;
        yA0             = 1-yB0;
        Evaporation_fxn2 = @(t, x) FuncEvaporation_OO(t, x,x30,yA0,yB0,dz,N,L,ColumnParam,CondParam,isotherm_parameters) ;

      
 %%% 4. SIMULATE CO2 WATER EVAPORATION STAGE [OPEN--OPEN]
        [t5, e] = ode15s(Evaporation_fxn2, [0 t_evap2], x30, opts4) ;

        % Correct Boundary Conditions 
        %%% Inlet
        e(:, 2*N+5)     = PressureVaporSat(Tf);
        idx             = find(e(:, 2*N+6) > e(:,2*N+5))    ; % Find values of P(2) > P(1)
        e(idx, 2*N+5)   = e(idx, 2*N+6)                     ;        
        e(:, 1)         = yA0                               ; 
        e(:, N+3)       = yB0                               ;
        e(:, 3*N+7)     = e(:, 3*N+8)                       ; % qCO2(1)   = qCO2(2)
        e(:, 4*N+9)     = e(:, 4*N+10)                      ; % qH20(1)    = qH20(2)                
        % e(:, 5*N+11)    = Tf; 
        
        MW_in           = MW_air + (MW_water-MW_air).*e(:, N+3);
        rho_g1     = e(:,2*N+5)/(R*Ta);  
        kinetic_term  = (1.75*(1-eb).*rho_g1.*MW_in)./(eb*dp) ;
        dsmt1   = viscous_term^2 - (4*kinetic_term.*(e(:,2*N+6)-e(:,2*N+5))*(2/dz)) ;
        u    = (-viscous_term + sqrt(dsmt1))./(2*kinetic_term) ;     
        
        rho_g           = e(:, 2*N+5)./(R*Tf)    ; % P(1)/(R*T(1))
        Constant2       = rho_g.*((eb.*u.*C_pg*dz)./(2*Kz)) ;
        e(:, 5*N+11)    = (e(:, 5*N+12) + (Tf.*Constant2(:)))./(1+Constant2(:));      % T(1) (Dankwerts BC)

        %%% Outlet
        e(:, N+2)       = e(:, N+1)                         ; % yCO2(N+2) = yCO2(N+1)        
        e(:, 2*N+4)     = e(:, 2*N+3)                       ; % yH20(N+2) = yCO2(N+1)
        e(:, 4*N+8)     = e(:, 4*N+7)                       ; % qCO2(N+2) = qCO2(N+1)
        e(:, 5*N+10)    = e(:, 5*N+9)                       ; % qH20(N+2)  = qH20(N+1)
        e(:, 6*N+12)    = e(:, 6*N+11)                      ;
        MW_out          = MW_air + (MW_water-MW_air).*e(:, 2*N+4)           ; % T(N+2) = T(N+1)        
        e(: ,3*N+6)     = (-(viscous_term*u_vp*dz2) + e(: ,3*N+5))./(1 + (dz/2/R./e(:, 5*N+11)).*(1.75*(1-eb).*MW_out/dp/eb).*u_vp^2);                                                   
        yCe             = 1-(e(:,1:N+2) + e(:,N+3:2*N+4))   ; % yN2(1:N+2)

        % Prepare initial conditions for CO2 Desorption Step (Closed-open)
        x40         = e(end, :)'    ; % Final state of previous step is the initial state for current step
       
        Desorption_fxn        = @(t, x) FuncDesorption(t, x, x40, dz,N,L,ColumnParam,CondParam,isotherm_parameters) ;

  %%% 5. SIMULATE CO2 REMOVAL STEP [CLOSED--OPEN]
        [t6, f] = ode15s(Desorption_fxn, [0 t_des], x40, opts5) ;

        % Correct the output (clean up results from simulation)
        %%% Inlet
        f(:, 2*N+5)     = f(:, 2*N+6)                   ; % P(1) = P(2)
        f(:, 1)         = f(:, 2)                       ; % yCO2(1) = yCO2(2)
        f(:, N+3)       = f(:, N+4)                     ; % yH20(1) = yH20(2)
        f(:, 3*N+7)     = f(:, 3*N+8)                   ; % qCO2(1) = qCO2(2)
        f(:, 4*N+9)     = f(:, 4*N+10)                  ; % qH20(1) = qH20(2)   
        f(:, 5*N+11)    = f(:, 5*N+12)                  ; % T(1) = T(2)
        
        %%% Outlet
        idx             = find(f(:, 3*N+5) < f(:,3*N+6))    ; % Find values of P(N+1) < P(N+2)
        f(idx, 3*N+6)   = f(idx, 3*N+5)    ;     
        f(:, N+2)       = f(:, N+1)                     ; % yCO2(N+2) = yCO2(N+1)       
        f(:, 2*N+4)     = f(:, 2*N+3)                   ; % yH20(N+2) = yH20(N+1)
        f(:, 4*N+8)     = f(:, 4*N+7)                   ; % qCO2(N+2) = qCO2(N+1)
        f(:, 5*N+10)    = f(:, 5*N+9)                   ; % qH20(N+2) = qH20(N+1)        
        f(:, 6*N+12)    = f(:, 6*N+11)                  ; % T(N+2) = T(N+1)
        yCf             = 1-(f(:,1:N+2) + f(:,N+3:2*N+4)) ;   % yN2(1:N+2)

        % Re-define mole fraction of water to original feed value. 
        RH              = RH2;
        PH2O_sat        = PressureVaporSat(Ta);
        yB0             = (RH*PH2O_sat)/(100*PH)      ;
        yA0             = 4e-4;


        % Prepare initial conditions for Pressurization Step
        x50         = f(end, :)'; % Final state of previous step is the initial state for current step
        
        % Final conditions of states at first step of the Cycle
        statesFC = f(end, [2:N+1, N+4:2*N+3, 2*N+6:3*N+5, 3*N+8:4*N+7, 4*N+10:5*N+9, 5*N+12:6*N+11]);        % the only unknowns we consider for CSS condition is yCO2,P,qCO2,T,Twall
       
        % Set initail conditions for the next step
        x0 = x50;

        % Check Cyclic Steady State
            dyCO2 = abs((statesFC(1:N) - statesIC(1:N))/statesIC(1:N));
            dyH2O = abs((statesFC(N+1:2*N) - statesIC(N+1:2*N))/statesIC(N+1:2*N));
            dPress = abs((statesFC(2*N+1:3*N) - statesIC(2*N+1:3*N))/statesIC(2*N+1:3*N));
            dqCO2 = abs((statesFC(3*N+1:4*N) - statesIC(3*N+1:4*N))/statesIC(3*N+1:4*N));
            dqH2O = abs((statesFC(4*N+1:5*N) - statesIC(4*N+1:5*N))/statesIC(4*N+1:5*N));
            dTemp = abs((statesFC(5*N+1:6*N) - statesIC(5*N+1:6*N))/statesIC(5*N+1:6*N));
            MassBalance = MassBalanceCheck_New(a,b,c,e,f,t1,t2,t3,t5,t6,N,L,ColumnParam,CondParam) ;
            mass_balance = MassBalance.CO2_cycle_balance ;            
            cycle_data(i,:) = [i, dyCO2, dyH2O, dPress, dqCO2, dqH2O, dTemp, abs(mass_balance-1)]; % Store cycle number, CSS, and mass balance error
            
      % Check for convergence based on tolerance for dyCO2, dyH2O, dPress, dqCO2, dqH2O, dTemp
        if  all(abs(dyCO2) < CSS_tolerance) && all(abs(dyH2O) < CSS_tolerance) && ...
            all(abs(dPress) < CSS_tolerance) && all(abs(dqCO2) < CSS_tolerance) && ...
            all(abs(dqH2O) < CSS_tolerance) && all(abs(dTemp) < CSS_tolerance) && ...            
            abs(mass_balance - 1) <= MBS_tolerance
            consecutive_cycles_met = consecutive_cycles_met + 1;
        else
            consecutive_cycles_met = 0; % Reset counter if tolerance is not met
        end

      % Check if the tolerance criteria were met for 5 consecutive cycles
        if consecutive_cycles_met >= 5
            fprintf('Breaking loop after %d cycles: All variables within tolerance for 5 consecutive cycles.\n', i);
            break_condition_met = true; % Mark condition as met
            break; % Exit the loop
        end
      
      % % Write cycle data to Excel after the loop ends
      %   cycle_data = cycle_data(1:i,:) ; % Update cycle data if only few cycles were used.
      %   filename = 'Cycle_Data.xlsx';           
      %   headers = {'Cycle Number', 'dyCO2', 'dyH2O', 'dPress', 'dqCO2', 'dqH2O', 'dTemp', 'MB_Criteria'};
      %   writecell(headers, filename, 'Sheet', 'Sheet1', 'Range', 'A1');
      %   writematrix(cycle_data, filename, 'Sheet', 'Sheet1', 'Range', 'A2');
end

%% Display Peformance Parameters
toc
% PERFORMANCE METRICS
  Objectives = MechanicalEnergy(a, b, c, e, f, t1, t2, t3, t5, t6, RH, N, L, ColumnParam, CondParam, Times);
  
  MECH_ENERGY       = Objectives(1) ; 
  CO2_PRODUCTIVITY  = Objectives(2);
  gH2Olost_gCO2cap  = Objectives(3);
  gH2Oused_gCO2cap   = Objectives(4);  
  CO2_Purity    = Objectives(5);
  

end