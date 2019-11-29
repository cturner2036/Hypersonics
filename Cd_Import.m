%Cd Importing Script
[spreadsheet,filepath] = uigetfile('Cd_ss.xlsx');
Aero_Data = strcat(filepath, spreadsheet);
global NFC
NFC = readtable(Aero_Data, 'Sheet', 'Aero Properties', 'Range', 'A1:I49', 'ReadVariableNames', 1, 'PreserveVariableNames', 1, 'UseExcel', false);
global AFC
AFC = readtable(Aero_Data, 'Sheet', 'Aero Properties', 'Range', 'A50:I98', 'ReadVariableNames', 1, 'PreserveVariableNames', 1, 'UseExcel', false);
disp("Reading in Values from Aero Coefficients Table for 5X GHV)