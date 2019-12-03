%MASTER SCRIPT%
%
%This is the scramjet master script. It calls a number of functions that
%runs a 1-D scram jet engine mission.
%%
%Inputs:
% -InletMap > obtain this value from 'InportMap' fn
% -DynamicPressure > Q value, can change with vehicle position, (maybe
%                    remove this input to a lower level fn.
% -FreestreamMach > Drop Mach value
% -AngleofAttack > Drop Angle of Attack
% -Dry_Weight > Dry weight of the vehicle
% -Fuel_Weight > starting weight of Fuel that vehicle can carry
% -SI_Flag > switches between SI and English units (0 or 1)
% [Need to have '3X GHV Inlet Map-1.xlsx' & 'Cd_ss.xlsx' spreadsheets in
% WD]
%Required Functions:
%[ScramjetEngine, getGHV, Station, Plug_Nozzle, Plug_Nozzle_Style2,
%Flow_Properties, DragCoeff, atmosphere, actualWeight]
%
%%
%Output:
%Video of the mission (maybe a multiplot video)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

%ImportData
GHVData = ImportMap(5);
Cd_Import;

%Initial Conditions
Q = 1500;   %initial Q
FreeMach = 4.0;
AngleofAttack = 0.0;
Dry_Weight = 5051.142*.453592;  %[kgs]
Fuel_Weight = 3141.649*.453592; %[kgs]
FP_Angle = 2;
SI_Flag = 0;

%Initiate Mission
EnginePerf = ScramjetEngine(GHVData,Q,FreeMach,AngleofAttack,Dry_Weight,Fuel_Weight,SI_Flag,AFC,NFC,FP_Angle);
