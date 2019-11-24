function GHV_Instance = getGHV(GHV_dataset,Mach,AoA,DynamicPressure,SI_flag)
%% GETGHV Function getGHV
% This functions takes the GHV dataset and reports a specific condition
% 
% Inputs include:
% 
% GHVdataset - this is the GHV inlet map in it's imported form already scaled 
% for global size
% 
% Mach - The flight Mach number condition
% 
% AoA - the flight angle of attack relative to the velocity
% 
% Dynamic Pressure - the function will scale from the dynamic pressure in 
% the data set to the flight dynamic pressure specified
% 
%         if SI = 0 then Dynamic Pressure is expected to be lb/ft²
% 
%         if SI = 1 then Dynamic Pressure is expected to be in Pascals or 
% N/m²
% 
% SI_flag - 0 specifies customery English units, 1 specifies SI units
xq = Mach;
yq = AoA;
x = GHV_dataset.x_Mach;
y = GHV_dataset.y_AoA;
%Conversion Constants
Pa_psi = 6894.8;
m_ft = 0.3048;
m_in = 0.0254;
degR_K = 1.8;
kg_lb = 0.4536;
kJ_kg_BTU_lbm = 2.326;
N_lb = 4.44822;
GHV_Instance.FlightMachNumber = interp2(x, y, GHV_dataset.FlightMachNumber, xq,yq);
GHV_Instance.AngleofAttack_deg = interp2(x, y, GHV_dataset.AngleofAttack_deg,xq,yq);
GHV_Instance.YaawAngle_deg = interp2(x, y, GHV_dataset.YawAngle_deg,xq,yq);
if (0 == SI_flag)
    GHV_Instance.DynamicPressure_psf = DynamicPressure;
    GHV_Instance.DynamicPressure_Pa = DynamicPressure*Pa_psi/144; % convert psf to psi to Pa
else
    GHV_Instance.DynamicPressure_pa = DynamicPressure;
    GHV_Instance.DynamicPressure_psf = DynamicPressure*144/Pa_psi; % convert Pa to psi to psf
end
Qscale = GHV_Instance.DynamicPressure_psf / interp2(x, y, GHV_dataset.DynamicPressure_psf,xq,yq);
[T0_degR, P0_psi, rho0_slugft3,hgeom_ft] = AtmosQM(GHV_Instance.DynamicPressure_psf,Mach);
T_scale = T0_degR/interp2(x, y, GHV_dataset.FreestreamTemperature_degR,xq,yq);
P_scale = P0_psi / interp2(x,y,GHV_dataset.FreestreamPressure_psia,xq,yq);
m_scale = P_scale / sqrt(T_scale);
%GHV_Instance.Altitude_ft = interp2(x, y, GHV_dataset.Altitude_ft,xq,yq);
GHV_Instance.Altitude_ft = hgeom_ft;
GHV_Instance.Altitude_m = GHV_Instance.Altitude_ft*m_ft; 
%GHV_Instance.FreestreamPressure_psia = interp2(x, y, GHV_dataset.FreestreamPressure_psia,xq,yq);
GHV_Instance.FreestreamPressure_psia = P0_psi;
GHV_Instance.FreestreamPressure_Pa=GHV_Instance.FreestreamPressure_psia*Pa_psi;
GHV_Instance.FreestreamTotalPressure_psia = interp2(x, y, GHV_dataset.FreestreamTotalPressure_psia,xq,yq)*P_scale;
GHV_Instance.FreestreamTotalPressure_Pa=GHV_Instance.FreestreamTotalPressure_psia*Pa_psi;
GHV_Instance.FreestreamTemperature_degR = interp2(x, y, GHV_dataset.FreestreamTemperature_degR,xq,yq)*T_scale;
GHV_Instance.FreestreamTemperature_K = GHV_Instance.FreestreamTemperature_degR / degR_K;
GHV_Instance.FreeStreamTotalTemperature_degR = interp2(x, y, GHV_dataset.FreestreamTotalTemperature_degR,xq,yq)*T_scale;
GHV_Instance.FreeStreamTotalTemperature_K = GHV_Instance.FreeStreamTotalTemperature_degR / degR_K;
GHV_Instance.FreestreamEquilibriumGama = interp2(x, y, GHV_dataset.FreestreamEquilibriumGamma,xq,yq);
GHV_Instance.FreestreamVelocity_fts = interp2(x, y, GHV_dataset.FreestreamVelocity_fts,xq,yq)*sqrt(T_scale);
GHV_Instance.FreestreamVelocity_ms = GHV_Instance.FreestreamVelocity_fts * m_ft;
GHV_Instance.PhysicalCaptureArea_in2 = interp2(x, y, GHV_dataset.PhysicalCaptureArea_in2,xq,yq);
GHV_Instance.PhysicalCaptureArea_m2 = GHV_Instance.PhysicalCaptureArea_in2 *(m_in^2);
GHV_Instance.EffectiveStreamtubeCapture_in2 = interp2(x, y, GHV_dataset.EffectiveStreamtubeCapture_in2,xq,yq);
GHV_Instance.EffectiveStreamtubeCapture_m2 = GHV_Instance.EffectiveStreamtubeCapture_in2 *(m_in^2);
GHV_Instance.InletMassCaptureRatio = interp2(x, y, GHV_dataset.InletMassCaptureRatio,xq,yq);
GHV_Instance.PhysicalContractionRatio = interp2(x, y, GHV_dataset.PhysicalContractionRatio,xq,yq);
GHV_Instance.ThroatArea_in2 = interp2(x, y, GHV_dataset.ThroatArea_in2,xq,yq);
GHV_Instance.ThroatArea_m2 = GHV_Instance.ThroatArea_in2 *(m_in^2);
GHV_Instance.EffectiveStreamtubeContractionRatio = interp2(x, y, GHV_dataset.EffectiveStreamtubeContractionRatio,xq,yq);
GHV_Instance.FreestreamMassFlowRate_lbms = interp2(x, y, GHV_dataset.FreestreamMassFlowRate_lbms,xq,yq)* m_scale;
GHV_Instance.FreestreamMassFlowRate_kgs = GHV_Instance.FreestreamMassFlowRate_lbms * kg_lb;
GHV_Instance.FreestreamMassFlowRate_lbms = GHV_Instance.FreestreamMassFlowRate_lbms*m_scale;
GHV_Instance.FreestreamMassFlowRate_kgs = GHV_Instance.FreestreamMassFlowRate_lbms*kg_lb;
GHV_Instance.MassFlowRate_lbms = interp2(x, y, GHV_dataset.CapturedMassFlowRate_lbms,xq,yq)*m_scale;
GHV_Instance.MassFlowRate_kgs = GHV_Instance.MassFlowRate_lbms * kg_lb;
GHV_Instance.BypassedAirflow_lbms = interp2(x, y, GHV_dataset.BypassedAirflow_lbms,xq,yq)*m_scale;
GHV_Instance.BypassedAirflow_kgs = GHV_Instance.BypassedAirflow_lbms * kg_lb;
GHV_Instance.KineticEnergyEfficiency = interp2(x, y, GHV_dataset.KineticEnergyEfficiency,xq,yq);
GHV_Instance.ThroatPressure_psia = interp2(x, y, GHV_dataset.ThroatPressure_psia,xq,yq)*P_scale;
GHV_Instance.ThroatPressure_Pa = GHV_Instance.ThroatPressure_psia * Pa_psi;
GHV_Instance.ThroatTotalPressure_psia = interp2(x, y, GHV_dataset.ThroatTotalPressure_psia,xq,yq)*P_scale;
GHV_Instance.ThroatTotalPressure_Pa = GHV_Instance.ThroatTotalPressure_psia * Pa_psi;
GHV_Instance.ThroatTemperature_degR = interp2(x, y, GHV_dataset.ThroatTemperature_degR,xq,yq)*T_scale;
GHV_Instance.ThroatTemperature_K = GHV_Instance.ThroatTemperature_degR /degR_K;
GHV_Instance.ThroatTotalTemperature_degR = interp2(x, y, GHV_dataset.ThroatTotalTemperature_degR,xq,yq)*T_scale;
GHV_Instance.ThroatTotalTemperature_K = GHV_Instance.ThroatTotalTemperature_degR / degR_K;
GHV_Instance.ThroatMachNumber = interp2(x, y, GHV_dataset.ThroatMachNumber,xq,yq);
GHV_Instance.ThroatGamma = interp2(x, y, GHV_dataset.ThroatGamma,xq,yq);
GHV_Instance.ThroatVelocity_fts = interp2(x, y, GHV_dataset.ThroatVelocity_fts,xq,yq)*sqrt(T_scale);
GHV_Instance.ThroatVelocity_ms = GHV_Instance.ThroatVelocity_fts *m_ft;
GHV_Instance.IsolatorHeatLoss_BTU_lbm = interp2(x, y, GHV_dataset.IsolatorHeatLoss_BTU_lbm,xq,yq); %likely scales with temperature, but not 1:1 assume small effect and ignore temperature scaling
GHV_Instance.IsolatorHeatLoss_kJ_kg = GHV_Instance.IsolatorHeatLoss_BTU_lbm * kJ_kg_BTU_lbm;
GHV_Instance.IsolatorExitArea_in2 = interp2(x, y, GHV_dataset.IsolatorExitArea_in2,xq,yq);
GHV_Instance.IsolatorExitArea_m2 = GHV_Instance.IsolatorExitArea_in2 * (m_in^2);
GHV_Instance.NormalShockPressureRatio = interp2(x, y, GHV_dataset.NormalShockPressureRatio,xq,yq);
GHV_Instance.IsolatorPressureRatio = interp2(x, y, GHV_dataset.IsolatorPressureRatio,xq,yq);
GHV_Instance.InletAxialForce_lbf = interp2(x, y, GHV_dataset.InletAxialForce_lbf,xq,yq)*P_scale; % m*V --> m_scale * sqrt(T_scale)=P_scale
GHV_Instance.InletAxialForce_N = GHV_Instance.InletAxialForce_lbf * N_lb;
GHV_Instance.InletAxialForceCoefficient = interp2(x, y, GHV_dataset.InletAxialForceCoefficient,xq,yq)*P_scale/Qscale;
% P_scale/Qscale should be 1, but differences in gamma assumption make this not so.
GHV_Instance.InletNormalForce_lbf = interp2(x, y, GHV_dataset.InletNormalForce_lbf,xq,yq)*P_scale;
GHV_Instance.InletNormalForce_N = GHV_Instance.InletNormalForce_lbf * N_lb;
GHV_Instance.InletNormalForceCoefficient = interp2(x, y, GHV_dataset.InletNormalForceCoefficient,xq,yq)*P_scale/Qscale;
end
