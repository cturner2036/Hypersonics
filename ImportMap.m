function GHV_dataset = ImportMap(GHV_size)
%% IMPORTMAP Import data from spreadsheet
% Function for importing data from the following spreadsheet:
%%
% 
%    Workbook: [path]\3X GHV Inlet Map-1.xlsx
%    Worksheet: 3X GHV Inlet
%
%% 
% Size is the relative scale for massflow.  Expected values range in integervalues 
% from 1-5 coresponding to the design massflow rate in 10-lb/sec increments
%% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 40);
% Specify sheet and range
opts.Sheet = "3X GHV Inlet";
opts.DataRange = "A4:AN51";
% Specify column names and types
opts.VariableNames = ["FlightMachNumber", "AngleofAttackdeg", "YawAngledeg", "DynamicPressurepsf", "Altitudeft", "FreestreamPressurepsia", "FreestreamTotalPressurepsia", "FreestreamTemperaturedegR", "FreestreamTotalTemperaturedegR", "FreestreamEquilibriumGamma", "FreestreamVelocityfts", "VarName12", "PhysicalCaptureAreain2", "EffectiveStreamtubeCapturein2", "InletMassCaptureRatio", "PhysicalContractionRatio", "ThroatAreain2", "EffectiveStreamtubeContractionRatio", "FreestreamMassFlowRatelbms", "CapturedMassFlowRatelbms", "BypassedAirflowlbms", "KineticEnergyEfficiency", "VarName23", "ThroatPressurepsia", "ThroatTotalPressurepsia", "ThroatTemperaturedegR", "ThroatTotalTemperaturedegR", "ThroatMachNumber", "ThroatGamma", "ThroatVelocityfts", "VarName31", "IsolatorHeatLossBTUlbm", "IsolatorExitAreain2", "NormalShockPressureRatio", "IsolatorPressureRatio", "VarName36", "InletAxialForcelbf", "InletAxialForceCoefficient", "InletNormalForcelbf", "InletNormalForceCoefficient"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "double", "double", "double", "double", "double", "double", "double", "string", "double", "double", "double", "double", "string", "double", "double", "double", "double"];
opts = setvaropts(opts, [12, 23, 31, 36], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [12, 23, 31, 36], "EmptyFieldRule", "auto");
% Import the data
[spreadsheet,filepath] = uigetfile('3X GHV Inlet Map-1.xlsx');
GHV_inlet_map = strcat(filepath, spreadsheet);
tbl = readtable(GHV_inlet_map, opts, "UseExcel", false);
%% Scale tbl to 5x vehicle
% input data is for 3X vehicle (30 lbs/sec)
% Target is GHV_size
GHV_dataset.size = GHV_size;
tmp_scale = GHV_dataset.size/3;
%tbl.FlightMachNumber;
%tbl.AngleofAttackdeg;
%tbl.YawAngledeg;
%tbl.DynamicPressurepsf;
%tbl.Altitudeft;
%tbl.FreestreamPressurepsia;
%tbl.FreestreamTotalPressurepsia;
%tbl.FreestreamTemperaturedegR;
%tbl.FreestreamTotalTemperaturedegR;
%tbl.FreestreamEquilibriumGamma;
%tbl.FreestreamVelocityfts;
tbl.PhysicalCaptureAreain2=tbl.PhysicalCaptureAreain2*tmp_scale;
tbl.EffectiveStreamtubeCapturein2=tbl.EffectiveStreamtubeCapturein2*tmp_scale;
%tbl.InletMassCaptureRatio;
%tbl.PhysicalContractionRatio;
tbl.ThroatAreain2=tbl.ThroatAreain2*tmp_scale;
%tbl.EffectiveStreamtubeContractionRatio;
tbl.FreestreamMassFlowRatelbms=tbl.FreestreamMassFlowRatelbms*tmp_scale;
tbl.CapturedMassFlowRatelbms=tbl.CapturedMassFlowRatelbms*tmp_scale;
tbl.BypassedAirflowlbms=tbl.BypassedAirflowlbms*tmp_scale;
%tbl.KineticEnergyEfficiency;
%tbl.ThroatPressurepsia;
%tbl.ThroatTotalPressurepsia;
%tbl.ThroatTemperaturedegR;
%tbl.ThroatTotalTemperaturedegR;
%tbl.ThroatMachNumber;
%tbl.ThroatGamma;
%tbl.ThroatVelocityfts;
%tbl.IsolatorHeatLossBTUlbm=tbl.IsolatorHeatLossBTUlbm;
tbl.IsolatorExitAreain2=tbl.IsolatorExitAreain2*tmp_scale;
%tbl.NormalShockPressureRatio;
%tbl.IsolatorPressureRatio;
tbl.InletAxialForcelbf=tbl.InletAxialForcelbf*tmp_scale;
tbl.InletAxialForceCoefficient=tbl.InletAxialForceCoefficient*tmp_scale;
tbl.InletNormalForcelbf=tbl.InletNormalForcelbf*tmp_scale;
%tbl.InletNormalForceCoefficient;
%% Convert to output type
GHV_dataset.FlightMachNumber = mytable(tbl.FlightMachNumber);
GHV_dataset.AngleofAttack_deg = mytable(tbl.AngleofAttackdeg);
GHV_dataset.YawAngle_deg = mytable(tbl.YawAngledeg);
GHV_dataset.DynamicPressure_psf = mytable(tbl.DynamicPressurepsf);
GHV_dataset.Altitude_ft = mytable(tbl.Altitudeft);
GHV_dataset.FreestreamPressure_psia = mytable(tbl.FreestreamPressurepsia);
GHV_dataset.FreestreamTotalPressure_psia = mytable(tbl.FreestreamTotalPressurepsia);
GHV_dataset.FreestreamTemperature_degR = mytable(tbl.FreestreamTemperaturedegR);
GHV_dataset.FreestreamTotalTemperature_degR = mytable(tbl.FreestreamTotalTemperaturedegR);
GHV_dataset.FreestreamEquilibriumGamma = mytable(tbl.FreestreamEquilibriumGamma);
GHV_dataset.FreestreamVelocity_fts = mytable(tbl.FreestreamVelocityfts);
GHV_dataset.PhysicalCaptureArea_in2 = mytable(tbl.PhysicalCaptureAreain2);
GHV_dataset.EffectiveStreamtubeCapture_in2 = mytable(tbl.EffectiveStreamtubeCapturein2);
GHV_dataset.InletMassCaptureRatio = mytable(tbl.InletMassCaptureRatio);
GHV_dataset.PhysicalContractionRatio = mytable(tbl.PhysicalContractionRatio);
GHV_dataset.ThroatArea_in2 = mytable(tbl.ThroatAreain2);
GHV_dataset.EffectiveStreamtubeContractionRatio = mytable(tbl.EffectiveStreamtubeContractionRatio);
GHV_dataset.FreestreamMassFlowRate_lbms = mytable(tbl.FreestreamMassFlowRatelbms);
GHV_dataset.CapturedMassFlowRate_lbms = mytable(tbl.CapturedMassFlowRatelbms);
GHV_dataset.BypassedAirflow_lbms = mytable(tbl.BypassedAirflowlbms);
GHV_dataset.KineticEnergyEfficiency = mytable(tbl.KineticEnergyEfficiency);
GHV_dataset.ThroatPressure_psia = mytable(tbl.ThroatPressurepsia);
GHV_dataset.ThroatTotalPressure_psia = mytable(tbl.ThroatTotalPressurepsia);
GHV_dataset.ThroatTemperature_degR = mytable(tbl.ThroatTemperaturedegR);
GHV_dataset.ThroatTotalTemperature_degR = mytable(tbl.ThroatTotalTemperaturedegR);
GHV_dataset.ThroatMachNumber = mytable(tbl.ThroatMachNumber);
GHV_dataset.ThroatGamma = mytable(tbl.ThroatGamma);
GHV_dataset.ThroatVelocity_fts = mytable(tbl.ThroatVelocityfts);
GHV_dataset.IsolatorHeatLoss_BTU_lbm = mytable(tbl.IsolatorHeatLossBTUlbm);
GHV_dataset.IsolatorExitArea_in2 = mytable(tbl.IsolatorExitAreain2);
GHV_dataset.NormalShockPressureRatio = mytable(tbl.NormalShockPressureRatio);
GHV_dataset.IsolatorPressureRatio = mytable(tbl.IsolatorPressureRatio);
GHV_dataset.InletAxialForce_lbf = mytable(tbl.InletAxialForcelbf);
GHV_dataset.InletAxialForceCoefficient = mytable(tbl.InletAxialForceCoefficient);
GHV_dataset.InletNormalForce_lbf = mytable(tbl.InletNormalForcelbf);
GHV_dataset.InletNormalForceCoefficient = mytable(tbl.InletNormalForceCoefficient);
GHV_dataset.x_Mach = GHV_dataset.FlightMachNumber(1,:);
GHV_dataset.y_AoA = GHV_dataset.AngleofAttack_deg(:,1);
GHV_dataset.SI_Flag = 0; % 0 = English units

% Add Constants
GHV_dataset.R_BTU_lbmolR = 1.986;
GHV_dataset.R_ftlbf_lbmolR = 1545;
GHV_dataset.MW_air = 28.965;
%Add Isolator Exit definitions
GHV_dataset.IsolatorExitPressure_psia = GHV_dataset.IsolatorPressureRatio .* GHV_dataset.ThroatPressure_psia;
GHV_dataset.IsolatorExitGamma = 1.4*ones(6,7);
GHV_dataset.IsolatorExitCp_BTU_lbmR = (GHV_dataset.R_BTU_lbmolR/GHV_dataset.MW_air)*GHV_dataset.IsolatorExitGamma./(GHV_dataset.IsolatorExitGamma-1);
GHV_dataset.IsolatorExitTotalTemperature_R = GHV_dataset.ThroatTotalTemperature_degR - GHV_dataset.IsolatorHeatLoss_BTU_lbm./GHV_dataset.IsolatorExitCp_BTU_lbmR;

% Find Effective Exit Area @Design
gc=32.174;
R2_ft_R = GHV_dataset.R_ftlbf_lbmolR / GHV_dataset.MW_air;
M2 = 1.46;
mygamma = GHV_dataset.IsolatorExitGamma(3,5);
Tt2 = GHV_dataset.IsolatorExitTotalTemperature_R(3,5);
T2_R = Tt2/(1+0.5*(mygamma-1)*M2^2);
P2_psia = GHV_dataset.IsolatorExitPressure_psia(3,5);
m_dot2 = GHV_dataset.CapturedMassFlowRate_lbms(3,5);
v2 = M2 * sqrt(mygamma*R2_ft_R*T2_R*gc);
A2_in2 = m_dot2 * R2_ft_R * T2_R / (P2_psia*v2);

%Fix effective exit area for all conditions
GHV_dataset.IsolatorEffectiveExitArea_in2 = A2_in2 * ones(6,7);

% Find isolator Exit Mach at all conditions - solve quadratic function of M^2
f1 = 1./(GHV_dataset.IsolatorExitGamma-1);
f2 = GHV_dataset.CapturedMassFlowRate_lbms.^2 .*GHV_dataset.IsolatorExitTotalTemperature_R.*R2_ft_R;
f3 = GHV_dataset.IsolatorEffectiveExitArea_in2.^2 .*GHV_dataset.IsolatorExitPressure_psia.^2 .*GHV_dataset.IsolatorExitGamma.*gc;
GHV_dataset.IsolatorExitMach = (-f1 + (f1.^2 + 2.*f1.*f2./f3).^0.5).^0.5;

% Find isolator Exit Temperature and total pressure
f1=(1+0.5*(GHV_dataset.IsolatorExitGamma-1).*GHV_dataset.IsolatorExitMach.^2);
GHV_dataset.IsolatorExitTemperature_R = GHV_dataset.IsolatorExitTotalTemperature_R./f1;
f2 = f1.^(GHV_dataset.IsolatorExitGamma./(GHV_dataset.IsolatorExitGamma-1));
GHV_dataset.IsolatorExitTotalPressure_psia = GHV_dataset.IsolatorExitPressure_psia.*f2;

end
%% Function to create 2D matrix of data
% Puts data in format needed to use interp2 function for querying vlaues using 
% Mach and AoA
% 
% use function interp2(x_Mach,y_AoA, Variable, xq,yq)
function output = mytable(input)
output = [input(1:6) input(8:13) input(15:20) input(22:27) input(29:34) input(36:41) input(43:48)];
end