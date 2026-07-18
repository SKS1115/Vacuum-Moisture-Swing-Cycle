% Thermal to Electrical Energy Conversion
% Input: Thermal energy in MJ/kgCO2
% Output: Electrical energy in kWh/tCO2

% Inputs
MJ_per_kgCO2 = [14.5, 16.4];   % Thermal energy values in MJ/kgCO2
efficiency = 0.8;         % Thermal-to-electrical conversion efficiency (e.g., 35%)

% Constants
MJ_to_kWh = 1 / 3.6;       % 1 kWh = 3.6 MJ
kg_to_tonne = 1000;        % 1 tonne = 1000 kg

% Thermal energy in kWh/tCO2
thermal_kWh_per_tCO2 = MJ_per_kgCO2 * MJ_to_kWh * kg_to_tonne;

% Electrical output in kWh/tCO2
electrical_kWh_per_tCO2 = thermal_kWh_per_tCO2 * efficiency;

% Display results
for i = 1:length(MJ_per_kgCO2)
    fprintf('%.2f MJ/kgCO2 (thermal) -> %.2f kWh/tCO2 (electrical) @ %.0f%% efficiency\n', ...
        MJ_per_kgCO2(i), electrical_kWh_per_tCO2(i), efficiency*100);
end
