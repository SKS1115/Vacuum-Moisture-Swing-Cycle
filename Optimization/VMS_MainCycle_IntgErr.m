function [MECH_ENERGY, CO2_PRODUCTIVITY, gH2Olost_gCO2cap, gH2Oused_gCO2cap, CO2_Purity] = VMS_MainCycle(N, RH, t_ad, t_vp, t_de, v_ads, v_vp)
% VMS_MainCycle (hardened against integration errors and array mismatches)
% - Parallel-safe: no globals, no persistent, no file I/O
% - Defensive ODE integration per step with early penalty returns
% - Robust CSS checks and size assertions before post-processing

tic;

%% --------- Problem constants / grids
L  = 0.01;
dz = L/N;
RH2 = RH;
dz2 = dz/2;

%% --------- Input parameters (pure-function style)
InputParameters     = ProcessInputParameters(N, RH, t_ad, t_vp, t_de, v_ads, v_vp);
ColumnParam         = InputParameters{1};
CondParam           = InputParameters{2};
isotherm_parameters = InputParameters{3};
Times               = InputParameters{4};

% Pull what we actually use
eb      = ColumnParam(5);
C_pg    = CondParam(3);
mu      = CondParam(7);
Dm      = CondParam(8);
Kz      = CondParam(10);
R       = CondParam(14);
dp      = CondParam(17);
PH      = CondParam(20);
PL      = CondParam(21);
v0      = CondParam(23);
Tf      = CondParam(24);
Ta      = CondParam(25);
yA0     = CondParam(26);
yB0     = CondParam(27);
u_vp    = CondParam(29);

MW_air   = 0.02896;
MW_water = 0.018015;
DL       = 0.7*Dm;
viscous_term = (150*(1-eb)^2*mu)/(eb^2*dp^2);

% Time slices for steps
t_pz    = Times(1);
t_ads   = Times(2);
t_evac  = Times(3);
t_evap1 = Times(4); %#ok<NASGU> % unused in provided code but kept for completeness
t_evap2 = Times(5);
t_des   = Times(6);

%% --------- Build ODE fun handles (pure)
Pressurization_fxn = @(t, x) FuncPressurization(t, x, dz, N, L, ColumnParam, CondParam, isotherm_parameters);
Adsorption_fxn     = @(t, x) FuncAdsorption(    t, x, dz, N, L, ColumnParam, CondParam, isotherm_parameters);
Evacuation_fxn     = @(t, x) FuncEvacuation(    t, x, dz, N, L, ColumnParam, CondParam, isotherm_parameters);
DesorptionCore_fxn = @(t, x) FuncDesorption(    t, x, dz, N, L, ColumnParam, CondParam, isotherm_parameters);
% Evaporation fun will be re-bound after step 3 because it depends on x30 and yA0/yB0

%% --------- Initial conditions (sanitized)
qstar = IsothermMS(yA0, yB0, PL, Ta, CondParam, isotherm_parameters);

x0                = zeros(7*N+14, 1);
x0(1:N+2)         = yA0;
x0(N+3:2*N+4)     = yB0;
x0(2*N+5:3*N+6)   = PL;
x0(3*N+7:4*N+8)   = qstar(1);
x0(4*N+9:5*N+10)  = qstar(2);
x0(5*N+11:6*N+12) = Ta;
x0(6*N+13:7*N+14) = Ta;

x0 = sanitizeIC(x0);

%% --------- ODE options (useful defaults for stiff systems)
baseOpts = odeset('RelTol', 1e-5, 'AbsTol', 1e-8, 'MaxStep', max(1e-3, 1e-4*max([t_pz, t_ads, t_evac, t_evap2, t_des])));
opts1 = baseOpts; opts2 = baseOpts; opts3 = baseOpts; opts4 = baseOpts; opts5 = baseOpts;

%% --------- Cycle loop / CSS
cycles = 30;
CSS_tolerance = 5e-3;
MBS_tolerance = 5e-3;
consecutive_cycles_met = 0;

% Pre-allocations (kept simple and robust)
cycle_data = zeros(cycles, 8);

% Initialize "previous cycle" subset of states for CSS test
statesIC = []; % will be set after first PZ

% Penalty to use on failure
PENALTY = struct( ...
    'MECH_ENERGY',        1e9, ...
    'CO2_PRODUCTIVITY',   0.0, ...
    'gH2Olost_gCO2cap',   1e9, ...
    'gH2Oused_gCO2cap',   1e9, ...
    'CO2_Purity',         0.0);

for i = 1:cycles

    %% 1) PRESSURIZATION [OPEN–CLOSED]
    [t1, a, ok1, msg1] = safeODE(Pressurization_fxn, [0, t_pz], x0, opts1);
    if ~ok1, [MECH_ENERGY, CO2_PRODUCTIVITY, gH2Olost_gCO2cap, gH2Oused_gCO2cap, CO2_Purity] = penaltyReturn(PENALTY, 'Pressurization', msg1); return; end
    a = applyBC_pressurization(a, N, eb, dp, R, Ta, yA0, yB0, C_pg, Kz, dz, dz2, viscous_term, DL, MW_air, MW_water);

    % Prepare initial IC for adsorption
    x10 = a(end, :)';
    x10(2*N+5:3*N+6) = PH; % set outlet pressure to PH
    x10 = sanitizeIC(x10);

    % CSS "previous" sample (subset of states at start of cycle)
    if isempty(statesIC)
        statesIC = a(1, [2:N+1, N+4:2*N+3, 2*N+6:3*N+5, 3*N+8:4*N+7, 4*N+10:5*N+9, 5*N+12:6*N+11]);
    end

    %% 2) ADSORPTION [OPEN–OPEN]
    [t2, b, ok2, msg2] = safeODE(Adsorption_fxn, [0, t_ads], x10, opts2);
    if ~ok2, [MECH_ENERGY, CO2_PRODUCTIVITY, gH2Olost_gCO2cap, gH2Oused_gCO2cap, CO2_Purity] = penaltyReturn(PENALTY, 'Adsorption', msg2); return; end
    b = applyBC_adsorption(b, N, eb, dp, R, Ta, yA0, yB0, PH, dz, dz2, viscous_term, DL, MW_air, MW_water);

    % Prepare IC for evacuation
    x20 = b(end, :)';
    x20 = sanitizeIC(x20);

    %% 3) EVACUATION [CLOSED–OPEN]
    [t3, c, ok3, msg3] = safeODE(Evacuation_fxn, [0, t_evac], x20, opts3);
    if ~ok3, [MECH_ENERGY, CO2_PRODUCTIVITY, gH2Olost_gCO2cap, gH2Oused_gCO2cap, CO2_Purity] = penaltyReturn(PENALTY, 'Evacuation', msg3); return; end
    c = applyBC_evacuation(c, N);

    % Prepare IC and inlet composition for evaporation
    x30 = c(end, :)'; x30 = sanitizeIC(x30);
    RH_local   = 100;
    Pvap_in    = PressureVaporSat(Tf);
    yB0_evap   = (RH_local*Pvap_in)/(100*Pvap_in);
    yA0_evap   = 1 - yB0_evap;

    Evaporation_fxn2 = @(t, x) FuncEvaporation_OO(t, x, x30, yA0_evap, yB0_evap, dz, N, L, ColumnParam, CondParam, isotherm_parameters);

    %% 4) EVAPORATION [OPEN–OPEN]
    [t5, e, ok4, msg4] = safeODE(Evaporation_fxn2, [0, t_evap2], x30, opts4);
    if ~ok4, [MECH_ENERGY, CO2_PRODUCTIVITY, gH2Olost_gCO2cap, gH2Oused_gCO2cap, CO2_Purity] = penaltyReturn(PENALTY, 'Evaporation', msg4); return; end
    e = applyBC_evaporation(e, N, eb, dp, R, Tf, Ta, C_pg, Kz, dz, dz2, viscous_term, u_vp, MW_air, MW_water);

    % Prepare IC for desorption
    x40 = e(end, :)'; x40 = sanitizeIC(x40);

    % Desorption uses x40 only to compute dependent terms in your original code; we keep the core fun pure
    Desorption_fxn = @(t, x) FuncDesorption(t, x, x40, dz, N, L, ColumnParam, CondParam, isotherm_parameters);

    %% 5) DESORPTION [CLOSED–OPEN]
    [t6, f, ok5, msg5] = safeODE(Desorption_fxn, [0, t_des], x40, opts5);
    if ~ok5, [MECH_ENERGY, CO2_PRODUCTIVITY, gH2Olost_gCO2cap, gH2Oused_gCO2cap, CO2_Purity] = penaltyReturn(PENALTY, 'Desorption', msg5); return; end
    f = applyBC_desorption(f, N);

    %% Reset feed (post-cycle)
    RH  = RH2; %#ok<NASGU>
    PH2O_sat = PressureVaporSat(Ta);
    yB0 = (RH2 * PH2O_sat) / (100 * PH);
    yA0 = 4e-4; %#ok<NASGU>

    %% Prepare next cycle IC
    x50 = f(end, :)';
    x0  = sanitizeIC(x50);

    %% CSS checks (robust to zero denominators)
    statesFC = f(end, [2:N+1, N+4:2*N+3, 2*N+6:3*N+5, 3*N+8:4*N+7, 4*N+10:5*N+9, 5*N+12:6*N+11]);

    epsDen = 1e-12;
    denom  = @(v) max(abs(v), epsDen);

    dyCO2  = abs((statesFC(1:N)               - statesIC(1:N)              ) ./ denom(statesIC(1:N)));
    dyH2O  = abs((statesFC(N+1:2*N)           - statesIC(N+1:2*N)          ) ./ denom(statesIC(N+1:2*N)));
    dPress = abs((statesFC(2*N+1:3*N)         - statesIC(2*N+1:3*N)        ) ./ denom(statesIC(2*N+1:3*N)));
    dqCO2  = abs((statesFC(3*N+1:4*N)         - statesIC(3*N+1:4*N)        ) ./ denom(statesIC(3*N+1:4*N)));
    dqH2O  = abs((statesFC(4*N+1:5*N)         - statesIC(4*N+1:5*N)        ) ./ denom(statesIC(4*N+1:5*N)));
    dTemp  = abs((statesFC(5*N+1:6*N)         - statesIC(5*N+1:6*N)        ) ./ denom(statesIC(5*N+1:6*N)));

    % Mass balance (pure)
    MassBalance  = MassBalanceCheck_New(a, b, c, e, f, t1, t2, t3, t5, t6, N, L, ColumnParam, CondParam);
    mass_balance = MassBalance.CO2_cycle_balance;

    cycle_data(i, :) = [i, mean(dyCO2), mean(dyH2O), mean(dPress), mean(dqCO2), mean(dqH2O), mean(dTemp), abs(mass_balance - 1)];

    if all(dyCO2 < CSS_tolerance) && all(dyH2O < CSS_tolerance) && ...
       all(dPress < CSS_tolerance) && all(dqCO2 < CSS_tolerance) && ...
       all(dqH2O < CSS_tolerance) && all(dTemp < CSS_tolerance) && ...
       abs(mass_balance - 1) <= MBS_tolerance
        consecutive_cycles_met = consecutive_cycles_met + 1;
    else
        consecutive_cycles_met = 0;
    end

    if consecutive_cycles_met >= 5
        % fprintf('Breaking loop after %d cycles: All variables within tolerance for 5 consecutive cycles.\n', i);
        break;
    end

    % Update statesIC for next CSS delta
    statesIC = statesFC;
end

%% --------- Post-processing: sanity assertions before MechanicalEnergy
% Ensure basic consistency so GA failures are caught early with clear errors
assert(size(a,2) == 7*N+14 && size(b,2) == 7*N+14 && size(c,2) == 7*N+14 && size(e,2) == 7*N+14 && size(f,2) == 7*N+14, ...
    'State matrices have unexpected number of columns. Check step functions output size.');

% Make time vectors column and strictly increasing
t1 = makeColVec(t1); t2 = makeColVec(t2); t3 = makeColVec(t3); t5 = makeColVec(t5); t6 = makeColVec(t6);

%% --------- PERFORMANCE METRICS (pure, no I/O)
Objectives = MechanicalEnergy(a, b, c, e, f, t1, t2, t3, t5, t6, RH, N, L, ColumnParam, CondParam, Times);

MECH_ENERGY       = Objectives(1);
CO2_PRODUCTIVITY  = Objectives(2);
gH2Olost_gCO2cap  = Objectives(3);
gH2Oused_gCO2cap  = Objectives(4);
CO2_Purity        = Objectives(5);

toc;

%% ====================== Nested helpers ======================

function x = sanitizeIC(x)
    % Replace NaN/Inf and enforce minimal positivity where needed
    x(~isfinite(x)) = 0;
    % Keep mole fractions in [0,1] if present (yA: 1..N+2, yB: N+3..2N+4)
    if numel(x) >= 2*N+4
        idxA = 1:N+2; idxB = N+3:2*N+4;
        x(idxA) = min(max(x(idxA), 0), 1);
        x(idxB) = min(max(x(idxB), 0), 1);
    end
    % Temperature floors
    if numel(x) >= 6*N+12
        idxT = 5*N+11:6*N+12;
        x(idxT) = max(x(idxT), 200); % K
    end
    % Pressure floors
    if numel(x) >= 3*N+6
        idxP = 2*N+5:3*N+6;
        x(idxP) = max(x(idxP), 1e-2); % Pa (avoid zero)
    end
end

function v = makeColVec(v)
    if isrow(v), v = v.'; end
end

function [t, y, ok, errmsg] = safeODE(fun, tspan, y0, opts)
    % Defensive integration wrapper
    try
        [t, y] = ode15s(fun, tspan, y0, opts);
        ok = true; errmsg = '';
        % Basic sanity check
        if isempty(t) || size(y,2) ~= numel(y0)
            ok = false; errmsg = 'Empty/size-mismatched ODE output.';
            t = tspan(:); y = repmat(y0(:).', numel(t), 1); % fallback shape
        end
        if any(~isfinite(y), 'all')
            ok = false; errmsg = 'Non-finite values in ODE states.';
        end
    catch ME
        ok = false; errmsg = ME.message;
        % Return minimally valid shapes so callers won't crash
        t = tspan(:);
        y = repmat(y0(:).', numel(t), 1);
    end
end

function [ME, CP, wLoss, wUse, pur] = penaltyReturn(P, stepName, why)
    % Optional: fprintf('Penalty on %s: %s\n', stepName, why);
    ME    = P.MECH_ENERGY;
    CP    = P.CO2_PRODUCTIVITY;
    wLoss = P.gH2Olost_gCO2cap;
    wUse  = P.gH2Oused_gCO2cap;
    pur   = P.CO2_Purity;
end

%% ---- Boundary condition fixups (pure algebra, no side effects)
function a = applyBC_pressurization(a, N, eb, dp, R, Ta, yA0, yB0, C_pg, Kz, dz, dz2, viscous_term, DL, MW_air, MW_water)
    % Inlet pressure monotone
    idx = a(:, 2*N+5) < a(:, 2*N+6);
    a(idx, 2*N+5) = a(idx, 2*N+6);
    MW_in = MW_air + (MW_water - MW_air).*a(:, N+3);
    rho_g1 = a(:, 2*N+5) / (R*Ta);
    kinetic_term = (1.75*(1-eb).*rho_g1.*MW_in) / (eb*dp);
    dsmt1 = viscous_term^2 - (4*kinetic_term .* (a(:, 2*N+6) - a(:, 2*N+5)) * (2/dz));
    u = max((-viscous_term + sqrt(max(dsmt1,0))) ./ (2*kinetic_term), 0);
    a(:, 1)     = (a(:,2)     + yA0 .* ((u.*dz) ./ (2*DL))) ./ (1 + ((u.*dz) ./ (2*DL)));
    a(:, N+3)   = (a(:,N+4)   + yB0 .* ((u.*dz) ./ (2*DL))) ./ (1 + ((u.*dz) ./ (2*DL)));
    rho_g = a(:, 2*N+5) / (R*Ta);
    Constant2 = rho_g .* ((eb.*u.*C_pg*dz) / (2*Kz));
    a(:, 5*N+11) = (a(:, 5*N+12) + (Ta .* Constant2(:))) ./ (1 + Constant2(:));
    % Zero-gradient at outlet
    a(:, N+2)     = a(:, N+1);
    a(:, 2*N+4)   = a(:, 2*N+3);
    a(:, 3*N+6)   = a(:, 3*N+5);
    a(:, 4*N+8)   = a(:, 4*N+7);
    a(:, 5*N+10)  = a(:, 5*N+9);
    a(:, 6*N+12)  = a(:, 6*N+11);
end

function b = applyBC_adsorption(b, N, eb, dp, R, Ta, yA0, yB0, PH, dz, dz2, viscous_term, DL, MW_air, MW_water)
    b(:, 5*N+11) = Ta;
    b(:, 1)      = (b(:,2)     + yA0*((dz) / (2*DL))*1*v0Corr()) ./ (1 + ((dz) / (2*DL))*1*v0Corr());
    b(:, N+3)    = (b(:,N+4)   + yB0*((dz) / (2*DL))*1*v0Corr()) ./ (1 + ((dz) / (2*DL))*1*v0Corr());
    % Pressure at inlet row using momentum correction
    MW_in = MW_air + (MW_water - MW_air) .* b(:, N+3);
    denom = 1 - (dz / (2*R) ./ b(:, 5*N+11)) .* (1.75*(1-eb).*MW_in/dp/eb) .* v0Corr().^2;
    b(:, 2*N+5) = ( (viscous_term*v0Corr()*dz2) + b(:, 2*N+6) ) ./ max(denom, 1e-12);
    % Zero-grad and PH at outlet
    b(:, 3*N+6) = PH;
    idx = b(:, 3*N+5) < PH;
    b(idx, 3*N+6) = b(idx, 3*N+5);
    b(:, N+2)     = b(:, N+1);
    b(:, 2*N+4)   = b(:, 2*N+3);
    b(:, 4*N+8)   = b(:, 4*N+7);
    b(:, 5*N+10)  = b(:, 5*N+9);
    b(:, 6*N+12)  = b(:, 6*N+11);
    function v = v0Corr()
        % If you embed v0 into CondParam, replace with that. Keeping v0 from base scope is ok.
        v = evalin('caller','v0'); %#ok<EVLDIR> % capture v0 from outer function.
    end
end

function c = applyBC_evacuation(c, N)
    % Inlet
    c(:, 2*N+5) = c(:, 2*N+6);
    c(:, 1)     = c(:, 2);
    c(:, N+3)   = c(:, N+4);
    c(:, 3*N+7) = c(:, 3*N+8);
    c(:, 4*N+9) = c(:, 4*N+10);
    c(:, 5*N+11)= c(:, 5*N+12);
    % Outlet
    idx = c(:, 3*N+6) > c(:, 3*N+5);
    c(idx, 3*N+6) = c(idx, 3*N+5);
    c(:, N+2)     = c(:, N+1);
    c(:, 2*N+4)   = c(:, 2*N+3);
    c(:, 4*N+8)   = c(:, 4*N+7);
    c(:, 5*N+10)  = c(:, 5*N+9);
    c(:, 6*N+12)  = c(:, 6*N+11);
end

function e = applyBC_evaporation(e, N, eb, dp, R, Tf, Ta, C_pg, Kz, dz, dz2, viscous_term, u_vp, MW_air, MW_water)
    e(:, 2*N+5) = PressureVaporSat(Tf);
    idx = e(:, 2*N+6) > e(:, 2*N+5);
    e(idx, 2*N+5) = e(idx, 2*N+6);
    % Inlet compositions fixed (from caller)
    % Keep given yA0_evap/yB0_evap by not overwriting here.
    e(:, 3*N+7) = e(:, 3*N+8);
    e(:, 4*N+9) = e(:, 4*N+10);
    MW_in = MW_air + (MW_water - MW_air).*e(:, N+3);
    rho_g1 = e(:, 2*N+5) / (R*Ta);
    kinetic_term = (1.75*(1-eb).*rho_g1.*MW_in)/(eb*dp);
    dsmt1 = viscous_term^2 - (4*kinetic_term .* (e(:,2*N+6) - e(:,2*N+5))*(2/dz));
    u = max((-viscous_term + sqrt(max(dsmt1,0))) ./ (2*kinetic_term), 0);
    rho_g = e(:, 2*N+5) / (R*Tf);
    Constant2 = rho_g .* ((eb.*u.*C_pg*dz)/(2*Kz));
    e(:, 5*N+11) = (e(:, 5*N+12) + (Tf.*Constant2(:))) ./ (1 + Constant2(:));
    % Outlet fixes
    e(:, N+2)     = e(:, N+1);
    e(:, 2*N+4)   = e(:, 2*N+3);
    e(:, 4*N+8)   = e(:, 4*N+7);
    e(:, 5*N+10)  = e(:, 5*N+9);
    e(:, 6*N+12)  = e(:, 6*N+11);
    MW_out = MW_air + (MW_water - MW_air) .* e(:, 2*N+4);
    denom = 1 + (dz/(2*R) ./ e(:, 5*N+11)) .* (1.75*(1-eb).*MW_out/dp/eb) .* u_vp.^2;
    e(:, 3*N+6) = (-(viscous_term*u_vp*dz2) + e(:, 3*N+5)) ./ max(denom, 1e-12);
end

function f = applyBC_desorption(f, N)
    % Inlet
    f(:, 2*N+5) = f(:, 2*N+6);
    f(:, 1)     = f(:, 2);
    f(:, N+3)   = f(:, N+4);
    f(:, 3*N+7) = f(:, 3*N+8);
    f(:, 4*N+9) = f(:, 4*N+10);
    f(:, 5*N+11)= f(:, 5*N+12);
    % Outlet
    idx = f(:, 3*N+5) < f(:, 3*N+6);
    f(idx, 3*N+6) = f(idx, 3*N+5);
    f(:, N+2)     = f(:, N+1);
    f(:, 2*N+4)   = f(:, 2*N+3);
    f(:, 4*N+8)   = f(:, 4*N+7);
    f(:, 5*N+10)  = f(:, 5*N+9);
    f(:, 6*N+12)  = f(:, 6*N+11);
end

end
