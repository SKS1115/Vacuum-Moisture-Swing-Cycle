
clear
clc 
tic 

% This script contains the main VMS cycle, computing the behavior of all 
% process steps—from pressurization through final desorption using ode15s. 
% The workflow is organized by defining the spatial domain, loading 
% operating-condition functions, retrieving required parameters, initializing
% the column state, and executing a step-wise loop across the VMS sequence. 
% The solver outputs matrices (a, b, c, d, e) for each step, representing 
% time (rows) and discretized volume L/dz (columns). State variables are 
% indexed as follows:
% 
% CO2 mole fraction  = (1:N+2) 
% H2O mole fraction  = (N+3:2*N+4) 
% Total Pressure     = (2*N+5:3*N+6) 
% CO2 loading        = (3*N+7:4*N+8) 
% H2O loading        = (4*N+9:5*N+10) 
% Total Temperature  = (5*N+11:6*N+12) 
% Wall Temperature   = (6*N+13:7*N+14) 
% N2 mole fraction   = 1 - ((1:N+2) + (N+3:2*N+4))


%``````````````````````````Domain Setup```````````````
% prompt1 = "column length (m) =";
% L = input(prompt1);             % Number of cells...number of rows (units: dimensionless)

dims = input('Enter [Inner Diameter (m), Outer Diameter (m), Bed Length (m)] = ');

if numel(dims) ~= 3
    error('Please enter exactly three values: [ID OD L]');
end

ID = dims(1);
OD = dims(2);
L  = dims(3);

prompt2 = "discretization size (-) =";
N = input(prompt2);             % Number of cells...number of rows (units: dimensionless)
prompt3 = "Relative Humidity (%) =" ;
RH = input(prompt3);
RH2 = RH;
prompt2 = "Number of Cycles (-) =" ;
cycles = input(prompt2);
cycle_data   = zeros(cycles,8) ; % This will hold the cycle number, CSS, and mass balance
break_condition_met = false    ; % Track whether the break condition is met
consecutive_cycles_met = 0     ;

% Ask user for key operating conditions all in ONE line
vals = input('Enter [vair(m/s) vvs(m/s) t_cs(s) t_vs(s) t_des(s)]: ');

% Basic validation
if numel(vals) ~= 5
    error('You must enter exactly 5 values: [vair vvs t_cs t_vs t_des]');
end

% Assign the input vals, their respective names
    vair  = vals(1);
    vvs   = vals(2);
    t_cs  = vals(3);
    t_vs  = vals(4);
    t_des = vals(5);

% Input parameters
    InputParameters     = ProcessInputParameters(N,RH,ID,OD,L,vair,vvs,t_cs,t_vs,t_des);   % Call input parameter function into this script of code
    ColumnParam         = InputParameters{1};               % Retrieve column parameters from input parameter function
    CondParam           = InputParameters{2};               % Retrieve condition parameters from input parameter function
    isotherm_parameters = InputParameters{3};               % Retrieve isotherm parameters from input parameter function
    Times               = InputParameters{4};               % Retrieve time parameters from input parameter function

% Retrieve Necessary parameters (pre-defined in "ProcessInputParameters") and re-assign respective names
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
    vair    = CondParam(23)  ;
    Tsat    = CondParam(24)  ;
    Tatm    = CondParam(25)  ;
    yA0     = CondParam(26)  ;
    yB0     = CondParam(27)  ;
    yC0     = CondParam(28)  ;
    vvs     = CondParam(29)  ;
    MW_air  = CondParam(34) ;
    MW_water= CondParam(35) ; 

% Define additional properties
    dz     = L/N            ; % column grid spacing
    dz2    = dz/2           ; % half column grid spacing
    DL     = 0.7*Dm         ; % axial dispersion 
    beta   = (150*(1-eb)^2*mu)/(eb^2*dp^2) ; % viscous term for ergun equation
    

% Create STEP function handles that embed all necessary parameters and state variables
    Pressurization_fxn    = @(t, x) FuncPressurization(t, x,dz,N,L,ColumnParam,CondParam,isotherm_parameters) ;
    CO2sorption_fxn       = @(t, x) FuncCO2sorption(t, x, dz,N,L,ColumnParam,CondParam,isotherm_parameters) ;
    Evacuation_fxn        = @(t, x) FuncEvacuation(t, x, dz,N,L,ColumnParam,CondParam,isotherm_parameters) ;
    Desorption_fxn        = @(t, x) FuncDesorption(t, x, dz,N,L,ColumnParam,CondParam,isotherm_parameters) ;

% Retrieve VMS Cycle Times
    t_pz    = Times(1);                       
    t_cs    = Times(2);                      
    t_ev    = Times(3); 
    t_vs    = Times(4); 
    t_des   = Times(5);  

% Column initialization (Initial Conditions (IC)) assuming low pressure
    % qstar               = IsothermMS(0.07, 0.93, PL, 275, CondParam, isotherm_parameters)  ; %(yCO2,yH20,Press, Temp, isotherm input)
    % x0                  = zeros(7*N+14,1) ;               % Initialize the size of the IC matrix
    % x0(1:N+2)           = 0.07             ;                 % Initial conditions for CO2 mole fraction
    % x0(N+3:2*N+4)       = 0.93           ;                 % Initial conditions for H20 mole fraction
    % x0(2*N+5:3*N+6)     = PL            ;                 % Initial conditions for total pressure
    % x0(3*N+7:4*N+8)     = qstar(1)      ;                 % Initial conditions for CO2 reaction rate
    % x0(4*N+9:5*N+10)    = qstar(2)      ;                 % Initial conditions for H20 reaction rate
    % x0(5*N+11:6*N+12)   = 275           ;                 % Initial conditions for total temperature
    % x0(6*N+13:7*N+14)   = 275           ;                 % Initial conditions for wall temperature

 % Column initialization (Initial Conditions (IC)) assuming high pressure
    qstar               = IsothermMS(yA0, yB0, PH, Tatm, CondParam, isotherm_parameters)  ; %(yCO2,yH20,Press, Temp, isotherm input)
    x0                  = zeros(7*N+14,1) ;                 % Initialize the size of the IC matrix
    x0(1:N+2)           = yA0             ;                 % Initial conditions for CO2 mole fraction
    x0(N+3:2*N+4)       = yB0             ;                 % Initial conditions for H20 mole fraction
    x0(2*N+5:3*N+6)     = PH              ;                 % Initial conditions for total pressure
    x0(3*N+7:4*N+8)     = qstar(1)        ;                 % Initial conditions for CO2 loading
    x0(4*N+9:5*N+10)    = qstar(2)        ;                 % Initial conditions for H20 loading
    x0(5*N+11:6*N+12)   = Tatm            ;                 % Initial conditions for total temperature
    x0(6*N+13:7*N+14)   = 298             ;                 % Initial conditions for wall temperature

%% Begin simulating VMS cycle.
%   Run the simulation until the change in the temperature, gas mole fraction
%   and CO2 molar loading is less than 0.5%. For the molar loading and the
%   mole fraction, the simulation also stops if the absolute change in the
%   state variable is less than the prescribed CSS criteria. 

    % Apply the options function to address tolerance issues in the ode solver
    opts1 = odeset( 'RelTol', 1e-3) ;
    opts2 = odeset( 'RelTol', 1e-3) ;
    opts3 = odeset( 'RelTol', 1e-3) ;
    opts4 = odeset( 'RelTol', 1e-3) ;
    opts5 = odeset( 'RelTol', 1e-3) ;


 for i = 1:cycles
     CSS_tolerance = 0.005  ; % Cyclic steady state criteria
     MBS_tolerance = 0.005  ; % Mass balance error criteria

     cycle_start = tic      ; % Track the time for each cycle

  %%% 1. SIMULATE PRESSURIZATION STEP [OPEN--CLOSED]
        [t1, a] = ode15s(Pressurization_fxn, [0 t_pz], x0, opts1) ;

        % Update the boundary conditions (clean up results from simulation)
        %%% Column Inlet
        idx             = find(a(:, 2*N+5) < a(:,2*N+5))                                    ; % Find values of P(1) < P(2)
        a(idx, 2*N+5)   = a(idx, 2*N+6)                                                     ; 
        MW_in           = MW_air + (MW_water-MW_air).*a(:, N+3)                             ; % molecular weight of species entering column
        rho_g           = a(:,2*N+5)/(R*Tatm)                                               ; % density in the gas phase
        alpha           = (1.75*(1-eb).*rho_g.*MW_in)./(eb*dp)                              ; % kinetic term for ergun equation
        disc            = beta^2 - (4*alpha.*(a(:,2*N+6)-a(:,2*N+5))*(2/dz))                ; % term used in ergun to compute velocity
        u               = (-beta + sqrt(disc))./(2*alpha)                                   ; % inlet velocity  
        a(:, 1)         = (a(:,2) + yA0.*((u.*dz)./(2*DL)))./ (1 + ((u.*dz)./(2*DL)))       ; % yCO2(1) = dankwerts BC
        a(:, N+3)       = (a(:,N+4) + yB0.*((u.*dz)./(2*DL)))./ (1 + ((u.*dz)./(2*DL)))     ; % yH2O(1) = dankwerts BC     
        convDiff        = rho_g.*((eb.*u.*C_pg*dz)./(2*Kz))                                 ; % term used in temperature BC
        a(:, 5*N+11)    = (a(:, 5*N+12) + (Tatm.*convDiff(:)))./(1+convDiff(:))             ; % T(1) = Dankwerts BC
        a(:, 3*N+7)     = a(:, 3*N+8)                                                       ; % qCO2(1) = qCO2(2)
        a(:, 4*N+9)     = a(:, 4*N+10)                                                      ; % qH20(1)  = H20(2)

        %%% Column Outlet
        a(:, N+2)       = a(:, N+1)                         ; % yCO2(N+2) = yCO2(N+1), imposing the zero gradient BC
        a(:, 2*N+4)     = a(:, 2*N+3)                       ; % yH20(N+2) = yH20(N+1), imposing the zero gradient BC
        a(: ,3*N+6)     = a(: ,3*N+5)                       ; % P(N+2) = P(N+1)
        a(:, 4*N+8)     = a(:, 4*N+7)                       ; % qCO2(N+2) = qCO2(N+1)
        a(:, 5*N+10)    = a(:, 5*N+9)                       ; % qH20(N+2)  = qH20(N+1)        
        a(:, 6*N+12)    = a(:, 6*N+11)                      ; % T(N+2) = T(N+1) 
        yCa = 1-(a(:,1:N+2) + a(:,N+3:2*N+4))               ; % yN2(1:N+2)
       
        % Prepare initial conditions for the CO2 sorption Step
        x10               = a(end, :)'  ; % Final state of pressurization becomes initial state for the next step
        x10(2*N+5:3*N+6)  = PH          ;

        % Save initial conditions of states at first step
        statesIC = a(1, [2:N+1, N+4:2*N+3, 2*N+6:3*N+5, 3*N+8:4*N+7, 4*N+10:5*N+9, 5*N+12:6*N+11]);  

  %%% 2. SIMULATE CO2 SORPTION STEP [OPEN--OPEN]
        [t2, b] = ode15s(CO2sorption_fxn, [0 t_cs], x10,opts2) ;
        
        % Update the boundary conditions (clean up results from simulation)
        %%% Column Inlet
        b(:, 1)         = (b(:,2) + yA0*((vair*dz)/(2*DL)))/ (1 + ((vair*dz)/(2*DL)))   ; % yCO2(1) = dankwerts BC
        b(:, N+3)       = (b(:,N+4) + yB0*((vair*dz)/(2*DL)))/ (1 + ((vair*dz)/(2*DL))) ; % yH20(1) = dankwerts BC
        rho_g           = PH/(R*Tatm)                                                   ; % gas density
        convDiff        = rho_g*((eb*vair*C_pg*dz)/(2*Kz))                              ; % term used in temp BC
        b(:, 5*N+11)    = (b(:, 5*N+12) + (Tatm.*convDiff(:)))./(1+convDiff(:))         ; % T(1) = Dankwerts BC     
        MW_in           = MW_air + (MW_water-MW_air)*b(:, N+3)                          ; % total molecular weight of species entering column
        b(:, 2*N+5)     = ((beta*vair*dz2) + b(:, 2*N+6))...
                          /(1-(dz/2/R/b(:, 5*N+11))*(1.75*(1-eb)*MW_in/dp/eb)*vair^2)   ; % P(1) = ergun
        b(:, 3*N+7)     = b(:, 3*N+8)                                                   ; % qCO2(1)   = qCO2(2)
        b(:, 4*N+9)     = b(:, 4*N+10)                                                  ; % qH20(1)    = qH20(2)

        %%% Column Outlet
        b(:, 3*N+6)     = PH                                ; % Fixed pressure PH   
        idx             = find(b(:, 3*N+5) < PH)            ; % find P(N+1) <PH = P(N+2)
        b(idx, 3*N+6)   = b(idx, 3*N+5)                     ; % P(N+2) = P(N+1)
        b(:, N+2)       = b(:, N+1)                         ; % yCO2(N+2) = yCO2(N+1)
        b(:, 2*N+4)     = b(:, 2*N+3)                       ; % yH20(N+2) = yCO2(N+1)     
        b(:, 4*N+8)     = b(:, 4*N+7)                       ; % qCO2(N+2) = qCO2(N+1)
        b(:, 5*N+10)    = b(:, 5*N+9)                       ; % qH20(N+2)  = qH20(N+1)        
        b(:, 6*N+12)    = b(:, 6*N+11)                      ; % T(N+2) = T(N+1)                   
        yCb             = 1-(b(:,1:N+2) + b(:,N+3:2*N+4))   ; % yN2(1:N+2)
              
        % Prepare initial conditions for N2 Evacuation Step
        x20         = b(end, :)'  ; % Final state of CO2 Sorption becomes initial state for the next step
        
     
  %%% 3. SIMULATE AIR EVACUATION STEP [CLOSED--OPEN]
        [t3, c] = ode15s(Evacuation_fxn, [0 t_ev], x20, opts3) ;
     
        % Update the boundary conditions (clean up results from simulation)
        %%% Column Inlet
        c(:, 2*N+5)     = c(:, 2*N+6)                       ; % P(1)   = P(2)
        c(:, 1)         = c(:, 2)                           ; % yCO2(1) = yCO2(2)
        c(:, N+3)       = c(:, N+4)                         ; % yH20(1) = yH20(2)
        c(:, 3*N+7)     = c(:, 3*N+8)                       ; % qCO2(1) = qCO2(2)
        c(:, 4*N+9)     = c(:, 4*N+10)                      ; % qH20(1) = qH20(2)
        c(:, 5*N+11)    = c(:, 5*N+12)                      ; % T(1) = T(2)
       
        %%% Column Outlet
        idx             = find(c(:, 3*N+6) > c(:,3*N+5))    ; % Find values of P(N+2) > P(N+1)
        c(idx, 3*N+6)   = c(idx, 3*N+5)                     ; % P(N+2)    = P(N+1) for P(N+1) values < P(N+2)
        c(:, N+2)       = c(:, N+1)                         ; % yCO2(N+2) = yCO2(N+1)
        c(:, 2*N+4)     = c(:, 2*N+3)                       ; % yH20(N+2) = yH20(N+1)       
        c(:, 4*N+8)     = c(:, 4*N+7)                       ; % qCO2(N+2) = qCO2(N+1)
        c(:, 5*N+10)    = c(:, 5*N+9)                       ; % qH20(N+2)  = qH20(N+1)
        c(:, 6*N+12)    = c(:, 6*N+11)                      ; % T(N+2) = T(N+1)
        yCc = 1-(c(:,1:N+2) + c(:,N+3:2*N+4))               ; % yN2(1:N+2)

        % Prepare initial conditions for Closed-Closed
        x30             = c(end, :)'     ; % Final state of air evacuation becomes initial state for the next step      

        % Initialize water vapor evaporation conditions
        RH           = 100                          ; % Water vapor relative humidity in step 4 (%)
        Psat_vs      = PressureVaporSat(Tsat)       ; % saturation pressure of water at Tsat
        yB0          = (RH*Psat_vs)/(100*Psat_vs)   ; % H2O inlet mole fraction in step 4
        yA0          = 1-yB0                        ; % CO2 inlet mole fraction in step 4

        % Create Vapor Stripping Function handle with embedded inputs
        VStripping_fxn = @(t, x) FuncVStripping(t, x,x30,yA0,yB0,dz,N,L,ColumnParam,CondParam,isotherm_parameters) ;

      
 %%% 4. SIMULATE CO2 WATER EVAPORATION STAGE [OPEN--OPEN]
        [t4, d] = ode15s(VStripping_fxn, [0 t_vs], x30, opts4) ;

        % Update the boundary conditions (clean up results from simulation)
        %%% Column Inlet
        d(:, 2*N+5)     = PressureVaporSat(Tsat)            ; % P(1)
        idx             = find(d(:, 2*N+6) > d(:,2*N+5))    ; % Find values of P(2) > P(1)
        d(idx, 2*N+5)   = d(idx, 2*N+6)                     ; % Update P(1)      
        d(:, 1)         = yA0                               ; % yCO2(1)
        d(:, N+3)       = yB0                               ; % yH2O(1) 
        d(:, 3*N+7)     = d(:, 3*N+8)                       ; % qCO2(1)  = qCO2(2)
        d(:, 4*N+9)     = d(:, 4*N+10)                      ; % qH20(1)  = qH20(2)                
        MW_in           = MW_air + (MW_water-MW_air).*d(:, N+3); %total molecular weight of species entering column
        rho_g           = d(:,2*N+5)/(R*Tsat)               ; % gas density
        alpha           = (1.75*(1-eb).*rho_g.*MW_in)./(eb*dp); % term used in ergun to compute velocity
        dpdz            = (d(:,2*N+6)-d(:,2*N+5))/(dz2)     ; % (P(2)-P(1))/dz/2
        disc            = beta^2 - (4*alpha.*dpdz)          ; % term used in ergun to compute velocity
        u               = (-beta + sqrt(disc))./(2*alpha)   ; % inlet velocity  
        convDiff      = rho_g.*((eb.*u.*C_pg*dz)./(2*Kz))   ;
        d(:, 5*N+11)    = (d(:, 5*N+12) + (Tsat.*convDiff(:)))./(1+convDiff(:));   % T(1) = Dankwerts BC

        %%% Column Outlet
        d(:, N+2)       = d(:, N+1)                         ; % yCO2(N+2) = yCO2(N+1)        
        d(:, 2*N+4)     = d(:, 2*N+3)                       ; % yH20(N+2) = yCO2(N+1)
        d(:, 4*N+8)     = d(:, 4*N+7)                       ; % qCO2(N+2) = qCO2(N+1)
        d(:, 5*N+10)    = d(:, 5*N+9)                       ; % qH20(N+2)  = qH20(N+1)
        d(:, 6*N+12)    = d(:, 6*N+11)                      ; % T(N+2) = T(N+1) 
        MW_out          = MW_air + (MW_water-MW_air).*d(:, 2*N+4); %total molecular weight of species leaving column      
        d(: ,3*N+6)     = (-(beta*vvs*dz2) + d(: ,3*N+5))...
                          ./(1 + (dz/2/R./d(:, 5*N+11)).*(1.75*(1-eb).*MW_out/dp/eb).*vvs^2); % P(N+2) = ergun                                             
        yCd             = 1-(d(:,1:N+2) + d(:,N+3:2*N+4))   ; % yN2(1:N+2)

        % Prepare initial conditions for CO2 Desorption Step (Closed-open)
        x40         = d(end, :)'    ; %Final state of vapor stripping becomes initial state for the next step
       
        Desorption_fxn        = @(t, x) FuncDesorption(t, x, x40, dz,N,L,ColumnParam,CondParam,isotherm_parameters) ;

  %%% 5. SIMULATE CO2 REMOVAL STEP [CLOSED--OPEN]
        [t5, e] = ode15s(Desorption_fxn, [0 t_des], x40, opts5) ;

        % Update the boundary conditions (clean up results from simulation)
        %%% Column Inlet
        e(:, 2*N+5)     = e(:, 2*N+6)                   ; % P(1) = P(2)
        e(:, 1)         = e(:, 2)                       ; % yCO2(1) = yCO2(2)
        e(:, N+3)       = e(:, N+4)                     ; % yH20(1) = yH20(2)
        e(:, 3*N+7)     = e(:, 3*N+8)                   ; % qCO2(1) = qCO2(2)
        e(:, 4*N+9)     = e(:, 4*N+10)                  ; % qH20(1) = qH20(2)   
        e(:, 5*N+11)    = e(:, 5*N+12)                  ; % T(1) = T(2)
        
        %%% Outlet
        idx             = find(e(:, 3*N+5) < e(:,3*N+6)); % Find values of P(N+1) < P(N+2)
        e(idx, 3*N+6)   = e(idx, 3*N+5)                 ;     
        e(:, N+2)       = e(:, N+1)                     ; % yCO2(N+2) = yCO2(N+1)       
        e(:, 2*N+4)     = e(:, 2*N+3)                   ; % yH20(N+2) = yH20(N+1)
        e(:, 4*N+8)     = e(:, 4*N+7)                   ; % qCO2(N+2) = qCO2(N+1)
        e(:, 5*N+10)    = e(:, 5*N+9)                   ; % qH20(N+2) = qH20(N+1)        
        e(:, 6*N+12)    = e(:, 6*N+11)                  ; % T(N+2) = T(N+1)
        yCf             = 1-(e(:,1:N+2) + e(:,N+3:2*N+4)) ; % yN2(1:N+2)

   % Re-define ambient air conditions for next cycle. 
        RH        = RH2                     ; % air RH
        PH2O_atm  = PressureVaporSat(Tatm)  ; % Partial pressure of air H2O
        yB0       = (RH*PH2O_atm)/(100*PH)  ; % H2O at air RH
        yA0       = 4e-4                    ; % CO2 at 400ppm

        % Prepare initial conditions for next Step
        x50      = e(end, :)'; % Final state of desorption becomes initial state for the next step
        x0       = x50      ;
        
        % Final conditions of states at the end of the Cycle
        statesFC = e(end, [2:N+1, N+4:2*N+3, 2*N+6:3*N+5, 3*N+8:4*N+7, 4*N+10:5*N+9, 5*N+12:6*N+11]);   
          
   % Cyclic Steady State Verification
        dyCO2  = abs((statesFC(1:N) - statesIC(1:N))/statesIC(1:N));
        dyH2O  = abs((statesFC(N+1:2*N) - statesIC(N+1:2*N))/statesIC(N+1:2*N));
        dPress = abs((statesFC(2*N+1:3*N) - statesIC(2*N+1:3*N))/statesIC(2*N+1:3*N));
        dqCO2  = abs((statesFC(3*N+1:4*N) - statesIC(3*N+1:4*N))/statesIC(3*N+1:4*N));
        dqH2O  = abs((statesFC(4*N+1:5*N) - statesIC(4*N+1:5*N))/statesIC(4*N+1:5*N));
        dTemp  = abs((statesFC(5*N+1:6*N) - statesIC(5*N+1:6*N))/statesIC(5*N+1:6*N));
        
        EnergyBalance   = EnergyBalanceFxn(a,b,c,d,e,t1,t2,t3,t4,t5,N,L,ColumnParam,CondParam,isotherm_parameters) ;
        MassBalance     = MassBalanceFxn(a,b,c,d,e,t1,t2,t3,t4,t5,N,L,ColumnParam,CondParam) ;
        
        mb_error = MassBalance.CO2_cycle_balance ; % extract mass balance error
        cycle_data(i,:) = round([i, dyCO2, dyH2O, dPress, dqCO2,dqH2O,dTemp, abs(mb_error-1)], 2, 'significant');  % Store cycle number, CSS, and mass balance error
        
    % Output cycle time and current cycle number
        cycle_time = toc(cycle_start);
        labels = {'Cycle', 'dyCO2', 'dyH2O', 'dPress', 'dqCO2', 'dqH2O', 'dTemp', 'MassBalanceErr'};
        cycletable = array2table(cycle_data(i, :), 'VariableNames', labels);      disp(cycletable);
    
    % Check for CSS convergence based on tolerance for dyCO2, dyH2O, dPress, dqCO2, dqH2O, dTemp
    if all(abs(dyCO2) < CSS_tolerance) && all(abs(dyH2O) < CSS_tolerance) && ...
       all(abs(dPress) < CSS_tolerance) && all(abs(dqCO2) < CSS_tolerance) && ...
       all(abs(dqH2O) < CSS_tolerance) && all(abs(dTemp) < CSS_tolerance) && ...
       abs(mb_error - 1) <= MBS_tolerance
        consecutive_cycles_met = consecutive_cycles_met + 1;
    else
        consecutive_cycles_met = 0;  % Reset counter if tolerance is not met
    end
    
    % Check if the CSS tolerance criteria were met for specified number of consecutive cycles
    if consecutive_cycles_met >= 5
        fprintf('\nBreaking loop after %d cycles: All variables within tolerance for 5 consecutive cycles.\n', i);
        break_condition_met = true;  % Mark condition as met
        break;  % Exit the loop
    end

    cycle_data = cycle_data(1:i,:);
    % save("TypeName.mat",'a','b','c','d','e','t1','t2','t3','t4','t5','L','N','-v7.3');
end

%% Display Peformance Parameters
last_cycle = i; 

disp('nCO2 Mass Balance Table'); disp(MassBalance.CO2_MassBalanceTable);
disp('nH2O Mass Balance Table'); disp(MassBalance.H2O_MassBalanceTable);

PerformanceMetrics = MechanicalEnergy(a, b, c, d, e, t1, t2, t3, t4, t5, RH, N, L, ColumnParam, CondParam, Times); 
toc
