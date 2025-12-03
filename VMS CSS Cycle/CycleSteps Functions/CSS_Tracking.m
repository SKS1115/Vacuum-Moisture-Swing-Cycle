% % clear all 
% % clc
% 
% % Load data from Excel
% filename = 'Cycle_Data.xlsx'; % Replace with your actual Excel file name
% data = readtable(filename);
% 
% % Extract data from the table
% cycle_number = data{:, 'CycleNumber'};  % Extract 'Cycle Number' column
% css = data{:, 'CSS'};  % Extract 'CSS' column
% mb_criteria = data{:, 'MB_Criteria'};  % Extract 'MB Criteria' column
% 
% % Create a figure for the plot
% figure;
% 
% % Plot CSS on the primary y-axis
% yyaxis left;  % This activates the left y-axis
% plot(cycle_number, css, '-o', 'LineWidth', 1.5, 'DisplayName', 'CSS');
% xlabel('Cycle Number');
% ylabel('CSS (<=0.005)');
% grid on;
% hold on;  % Hold the current plot to overlay the next plot
% 
% % Plot MB Criteria on the secondary y-axis
% yyaxis right;  % This activates the right y-axis
% plot(cycle_number, mb_criteria, '-s', 'LineWidth', 1.5, 'DisplayName', 'MB Criteria');
% ylabel('MB Criteria (<= 0.005)');
% 
%  % xlim([2 cycle_number(end)])
% % Add title and legend
% title('Cycle Data Analysis');
% % legend('show', 'Location', 'best');


%% 
% Read data from Excel sheet 'cycle_data.xlsx'
filename = 'Cycle_Data.xlsx'; % Replace with the actual path if needed
data = readtable(filename);  % Adjust the sheet name if necessary

% Extracting the relevant columns from the table
cycle_number = data{:, 'CycleNumber'};  % Replace 'Cycle_Number' with the actual column name
dyCO2 = data{:, 'dyCO2'};  % Replace 'dyCO2' with the actual column name
dyH2O = data{:, 'dyH2O'};  % Replace 'dyH2O' with the actual column name
dPress = data{:, 'dPress'}; % Replace 'dPress' with the actual column name
dqCO2 = data{:, 'dqCO2'}; % Replace 'dqCO2' with the actual column name
dqH2O = data{:, 'dqH2O'};  % Replace 'dqH2O' with the actual column name
dTemp = data{:, 'dTemp'};  % Replace 'dTemp' with the actual column name
MB_Criteria = data{:, 'MB_Criteria'};  % Replace 'MB_Criteria' with the actual column name

% Plotting each variable with a different color
figure;
hold on; % Keep all plots on the same figure

plot(cycle_number, dyCO2, '-o', 'DisplayName', 'yCO2', 'LineWidth', 2);
plot(cycle_number, dyH2O, '-s', 'DisplayName', 'yH2O', 'LineWidth', 2);
plot(cycle_number, dPress, '-^', 'DisplayName', 'Press', 'LineWidth', 2);
plot(cycle_number, dqCO2, '-d', 'DisplayName', 'qCO2', 'LineWidth', 2);
plot(cycle_number, dqH2O, '-x', 'DisplayName', 'qH2O', 'LineWidth', 2);
plot(cycle_number, dTemp, '-+', 'DisplayName', 'Temp', 'LineWidth', 2);
plot(cycle_number, MB_Criteria, '-*', 'DisplayName', 'MB Criteria', 'LineWidth', 2);

% Customize the plot
xlabel('Cycle Number');
ylabel('Tolerance');
title('State Variable CCS Tolerance<=0.005');
legend('show'); % Display the legend
grid on;
hold off; % Release the plot hold
