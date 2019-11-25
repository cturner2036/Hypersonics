myGHV = ImportMap(5);
x = myGHV.x_Mach;
y = myGHV.y_AoA;
createfigure(x,y,myGHV.YawAngle_deg, 'Yaw Angle', 'Angle (\circ)');
createfigure(x,y,myGHV.DynamicPressure_psf,'Dynamic Pressure','Pressure (psf)');
createfigure(x,y,myGHV.Altitude_ft, 'Altitude', 'Altitude (ft)');
createfigure(x,y,myGHV.FreestreamPressure_psia, 'Freestream Pressure', 'Pressure (psia)');
createfigure(x,y,myGHV.FreestreamTotalPressure_psia, 'Freestream Total Pressure', 'Pressure (psia)');
createfigure(x,y,myGHV.FreestreamTemperature_degR, 'Freestream Temperature', 'Temperature (\circR)');
createfigure(x,y,myGHV.FreestreamTotalTemperature_degR,'Freestream Total Temp.', 'Temperature (\circR)');
createfigure(x,y,myGHV.FreestreamEquilibriumGamma, 'Freestream Eq. \gamma','\gamma');
createfigure(x,y,myGHV.FreestreamVelocity_fts,'Freestream Velocity','Velocity (ft/sec)');
createfigure(x,y,myGHV.PhysicalCaptureArea_in2,'Physical Capture Area', 'Area (in^2)');
createfigure(x,y,myGHV.EffectiveStreamtubeCapture_in2, 'Effective Streamtube Capture', 'Area (in^2)');
createfigure(x,y,myGHV.InletMassCaptureRatio,'Inlet Mass Capture Ratio', 'Ratio (in^2/in^2)');
createfigure(x,y,myGHV.PhysicalContractionRatio,'Physical Contraction Ratio','Ratio (in^2/in^2)');
createfigure(x,y,myGHV.ThroatArea_in2, 'Isolator Throat Area', 'Area (in^2)');
createfigure(x,y,myGHV.EffectiveStreamtubeContractionRatio,'Eff. Streamtube Contraction Ratio', 'Ratio (in^2/in^2)');
createfigure(x,y,myGHV.FreestreamMassFlowRate_lbms,'Freestream Mass Flow Rate','Mass Flow (lb/sec)');
createfigure(x,y,myGHV.CapturedMassFlowRate_lbms,'Captured Mass Flow Rate','Mass Flow (lb/sec)');
createfigure(x,y,myGHV.BypassedAirflow_lbms,'Bypassed Air Flow', 'Mass Flow (lb/sec)');
createfigure(x,y,myGHV.KineticEnergyEfficiency, 'Kinetic Energy Efficiency', 'Efficiency');
createfigure(x,y,myGHV.ThroatPressure_psia,'Isolator Throat Pressure', 'Pressure (psia)');
createfigure(x,y,myGHV.ThroatTotalPressure_psia,'Isolator Total Throat Pressure', 'Pressure (psia)');
createfigure(x,y,myGHV.ThroatTemperature_degR, 'Isolator Throat Temp.', 'Temperature (\circR)');
createfigure(x,y,myGHV.ThroatTotalTemperature_degR, 'Isolator Total Throat Temp.', 'Temperature (\circR)');
createfigure(x,y,myGHV.ThroatMachNumber, 'Isolator Throat Mach Number', 'Throat Mach');
createfigure(x,y,myGHV.ThroatGamma, 'Isolator Throat \gamma', '\gamma');
createfigure(x,y,myGHV.ThroatVelocity_fts, 'Isolator Throat Velocity', 'Velocity (ft/sec)');
createfigure(x,y,myGHV.IsolatorHeatLoss_BTU_lbm, 'Isolator Heat Loss', 'Heat Loss (BTU/lb_m');
createfigure(x,y,myGHV.IsolatorExitArea_in2, 'Isolator Exit Area', 'Area (in^2)');
createfigure(x,y,myGHV.NormalShockPressureRatio,'Normal Shock Pressure Ratio', 'Ratio (psia/psia)');
createfigure(x,y,myGHV.IsolatorPressureRatio, 'Isolator Pressure Ratio', 'Ratio (psia/psia)');
createfigure(x,y,myGHV.InletAxialForce_lbf, 'Inlet Axial Force', 'Force (lb_f)');
createfigure(x,y,myGHV.InletAxialForceCoefficient, 'Inlet Axial Force Coefficient', 'Coefficient');
createfigure(x,y,myGHV.InletNormalForce_lbf, 'Normal Force', 'Force (lb_f)');
createfigure(x,y,myGHV.InletNormalForceCoefficient,'Normal Force Coefficient','Coefficient');

%%
function createfigure(xdata, ydata, zdata, myTitle, ztitle)
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.12 0.54 0.36 0.39]);
hold(axes1,'on');

% Create surf
mesh(xdata,ydata,zdata,'Parent',axes1);
xlabel('Mach');
ylabel('AoA (\circ)')
zlabel(ztitle)
title(strcat("Raw ",myTitle));

view(axes1,[-37.5 30]);
grid(axes1,'on');

% Create axes
axes2 = axes('Parent',figure1,...
    'Position',[0.12 0.06 0.36 0.39]);
hold(axes2,'on');
x = 4:0.1:7;
y = -4:0.2:6;
x2D=repelem(x,51,1);
y2D = repelem(y(:),1,31);
myLinearMesh = interp2(xdata,ydata,zdata,x2D,y2D,'linear');

% Create surf
mesh(x,y,myLinearMesh,'Parent',axes2);
xlabel('Mach');
ylabel('AoA (\circ)');
zlabel(ztitle);
title(strcat("Interpolated ", myTitle));

view(axes2,[-37.5 30]);
grid(axes2,'on');
% Create axes
axes3 = axes('Parent',figure1,...
    'Position',[0.60 0.54 0.36 0.39]);
hold(axes3,'on');
x3 = 4:0.1:8;
y3 = -4:0.2:6;
x2D=repelem(x3,51,1);
y2D = repelem(y3(:),1,41);
mySplineMesh = interp2(x,y,myLinearMesh,x2D,y2D,'spline');

% Create surf
mesh(x3,y3,mySplineMesh,'Parent',axes3);
xlabel('Mach');
ylabel('AoA (\circ)');
zlabel(ztitle);
title(strcat("Extrapolated ",myTitle));

view(axes3,[-37.5 30]);
grid(axes3,'on');

end