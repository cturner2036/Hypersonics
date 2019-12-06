function [qr, qw, Qr, Qw] = FullHeat(P0, T, M, Tw)
%{
ASEN 537 Final Project
Heat Transfer Script

Purpose: Compute the heat transfer required to keep the walls of the GHV w/
EERJ at a given temperature

Created: 11/24/19 19:23
%}


% assume max wall temp is 1000K
% P0 = 5.53e3;
R = 287;
a = 0;
% M0=4;M1=2;M2=1.7;M4=1;M9=8;
% T0=216.65;T1=250;T2=300;T4=2200;T9=1000;

Rho0 = P0/R/T(1);

Tw_Tip=Tw(1);
Tw_Inlet=Tw(2);
Tw_Isolator=Tw(3);
Tw_Combustor=Tw(4);
Tw_Ramp=Tw(5);

T0 = T(1); T1 = T(2); T2 = T(3); T4 = T(4); T9 = T(5);
M0 = M(1); M1 = M(2); M2 = M(3); M4 = M(4); M9 = M(5);


% define loaction vector x
n1=10;n2=100;n3=100;n4=100;n5=100;  

%set locations of stations
x_Inlet=1/100; x_Isolator = 1; x_Combustor = 2; x_Ramp=3; x_End = 5;%vehicle length in [m]

%% Section: Tip
%Fay and Riddell method for leading edge heat transfer
x1 = linspace(0,x_Inlet,n1); delX = (x_Inlet-0)/n1;
Width = ones(1,n1)*1; M = ones(1,n1)*M0;
for i = 1:length(x1)
    [q_Tip_r(i), q_Tip_w(i)] = LeadingEdgeHeat(M(i), a,T0,P0,Rho0, Tw_Tip);
    Q_Tip_r(i) = Width(i)*delX*sum(q_Tip_r);
    Q_Tip_w(i) = Width(i)*delX*sum(q_Tip_w);
end

%% Section: Inlet
% Eckhert Enthalpy method for for duct flow
T = linspace(T0,T1,n2);
x2 = linspace(x_Inlet, x_Isolator, n2); delX = (x_Isolator - x_Inlet)/n2;
Width = ones(1,n2); M = linspace(M0,M1,n2);
for i = 1:length(x2)
    [q_Inlet_r(i), q_Inlet_w(i) ]= DuctHeat(M(i),T(i),x2(i),Tw_Inlet);
    Q_Inlet_r(i) = Width(i)*delX*sum(q_Inlet_r);
    Q_Inlet_w(i) = Width(i)*delX*sum(q_Inlet_w);
end
% Christophel reference for shock induced b.l. separation hot-spots
%hmax/hinf = 1.1*(Pmax/Pinf)^0.7

%% Section: Isolator
T = linspace(T1,T2,n3);
x3 = linspace(x_Isolator, x_Combustor, n3); delX = (x_Combustor - x_Isolator)/n3;
Width = ones(1,n3); M = linspace(M1,M2,n3);
for i = 1:length(x3)
    [q_Isolator_r(i), q_Isolator_w(i)] = DuctHeat(M(i),T(i),x3(i),Tw_Isolator);
    Q_Isolator_r(i) = Width(i)*delX*sum(q_Isolator_r);
    Q_Isolator_w(i) = Width(i)*delX*sum(q_Isolator_w);
end

%% Section: Combustor
T = linspace(T2,T4,n4);
x4 = linspace(x_Combustor, x_Ramp, n4); delX = (x_Ramp - x_Combustor)/n4;
Width = ones(1,n4); M = linspace(M2,M4,n4);
for i = 1:length(x4)
    [q_Combustor_r(i), q_Combustor_w(i)] = DuctHeat(M(i),T(i),x4(i),Tw_Combustor);
    Q_Combustor_r(i) = Width(i)*delX*sum(q_Combustor_r);
    Q_Combustor_w(i) = Width(i)*delX*sum(q_Combustor_w);
    
end

%% Section: Exit Ramp
T = linspace(T4,T9,n5);
x5 = linspace(x_Ramp, x_End, n5); delX = (x_End - x_Ramp)/n5;
Width = ones(1,n5); M = linspace(M4,M9,n5);
for i = 1:length(x5)
    [q_Ramp_r(i), q_Ramp_w(i)] = FlatPlateHeat(M(i),a,T(i),x5(i),Tw_Ramp);
    Q_Ramp_r(i) = Width(i)*delX*sum(q_Ramp_r);
    Q_Ramp_w(i) = Width(i)*delX*sum(q_Ramp_w);    
end

% Rho-Mu High-Speed method
x = [x1 x2 x3 x4 x5]; 
qr = [q_Tip_r q_Inlet_r q_Isolator_r q_Combustor_r q_Ramp_r]; 
Qr = [Q_Tip_r Q_Inlet_r+Q_Tip_r(end) Q_Isolator_r+Q_Inlet_r(end)+Q_Tip_r(end)  Q_Combustor_r+Q_Isolator_r(end)+Q_Inlet_r(end)+Q_Tip_r(end) Q_Ramp_r+Q_Combustor_r(end)+Q_Isolator_r(end)+Q_Inlet_r(end)+Q_Tip_r(end)];
qw = [q_Tip_w q_Inlet_w q_Isolator_w q_Combustor_w q_Ramp_w]; 
Qw = [Q_Tip_w Q_Inlet_w+Q_Tip_w(end) Q_Isolator_w+Q_Inlet_w(end)+Q_Tip_w(end)  Q_Combustor_w+Q_Isolator_w(end)+Q_Inlet_w(end)+Q_Tip_w(end) Q_Ramp_w+Q_Combustor_w(end)+Q_Isolator_w(end)+Q_Inlet_w(end)+Q_Tip_w(end)];
end
% figure 
% hold on
% plot(x,Qr-Qw,x,Qr,x,Qw)
% figure
% plot(x,qr-qw,x,qr,x,qw)
=======
function [qr, qw, Qr, Qw] = FullHeat(P0, T, M, Tw)
%{
ASEN 537 Final Project
Heat Transfer Script

Purpose: Compute the heat transfer required to keep the walls of the GHV w/
EERJ at a given temperature

Created: 11/24/19 19:23
%}


% assume max wall temp is 1000K
% P0 = 5.53e3;
R = 287;
a = 0;
% M0=4;M1=2;M2=1.7;M4=1;M9=8;
% T0=216.65;T1=250;T2=300;T4=2200;T9=1000;

Rho0 = P0/R/T(1);

Tw_Tip=Tw(1);
Tw_Inlet=Tw(2);
Tw_Isolator=Tw(3);
Tw_Combustor=Tw(4);
Tw_Ramp=Tw(5);

T0 = T(1); T1 = T(2); T2 = T(3); T4 = T(4); T9 = T(5);
M0 = M(1); M1 = M(2); M2 = M(3); M4 = M(4); M9 = M(5);


% define loaction vector x
n1=10;n2=100;n3=100;n4=100;n5=100;  

%set locations of stations
x_Inlet=1/100; x_Isolator = 1; x_Combustor = 2; x_Ramp=3; x_End = 5;%vehicle length in [m]

%% Section: Tip
%Fay and Riddell method for leading edge heat transfer
x1 = linspace(0,x_Inlet,n1); delX = (x_Inlet-0)/n1;
Width = ones(1,n1)*1; M = ones(1,n1)*M0;
for i = 1:length(x1)
    [q_Tip_r(i), q_Tip_w(i)] = LeadingEdgeHeat(M(i), a,T0,P0,Rho0, Tw_Tip);
    Q_Tip_r(i) = Width(i)*delX*sum(q_Tip_r);
    Q_Tip_w(i) = Width(i)*delX*sum(q_Tip_w);
end

%% Section: Inlet
% Eckhert Enthalpy method for for duct flow
T = linspace(T0,T1,n2);
x2 = linspace(x_Inlet, x_Isolator, n2); delX = (x_Isolator - x_Inlet)/n2;
Width = ones(1,n2); M = linspace(M0,M1,n2);
for i = 1:length(x2)
    [q_Inlet_r(i), q_Inlet_w(i) ]= DuctHeat(M(i),T(i),x2(i),Tw_Inlet);
    Q_Inlet_r(i) = Width(i)*delX*sum(q_Inlet_r);
    Q_Inlet_w(i) = Width(i)*delX*sum(q_Inlet_w);
end
% Christophel reference for shock induced b.l. separation hot-spots
%hmax/hinf = 1.1*(Pmax/Pinf)^0.7

%% Section: Isolator
T = linspace(T1,T2,n3);
x3 = linspace(x_Isolator, x_Combustor, n3); delX = (x_Combustor - x_Isolator)/n3;
Width = ones(1,n3); M = linspace(M1,M2,n3);
for i = 1:length(x3)
    [q_Isolator_r(i), q_Isolator_w(i)] = DuctHeat(M(i),T(i),x3(i),Tw_Isolator);
    Q_Isolator_r(i) = Width(i)*delX*sum(q_Isolator_r);
    Q_Isolator_w(i) = Width(i)*delX*sum(q_Isolator_w);
end

%% Section: Combustor
T = linspace(T2,T4,n4);
x4 = linspace(x_Combustor, x_Ramp, n4); delX = (x_Ramp - x_Combustor)/n4;
Width = ones(1,n4); M = linspace(M2,M4,n4);
for i = 1:length(x4)
    [q_Combustor_r(i), q_Combustor_w(i)] = DuctHeat(M(i),T(i),x4(i),Tw_Combustor);
    Q_Combustor_r(i) = Width(i)*delX*sum(q_Combustor_r);
    Q_Combustor_w(i) = Width(i)*delX*sum(q_Combustor_w);
    
end

%% Section: Exit Ramp
T = linspace(T4,T9,n5);
x5 = linspace(x_Ramp, x_End, n5); delX = (x_End - x_Ramp)/n5;
Width = ones(1,n5); M = linspace(M4,M9,n5);
for i = 1:length(x5)
    [q_Ramp_r(i), q_Ramp_w(i)] = FlatPlateHeat(M(i),a,T(i),x5(i),Tw_Ramp);
    Q_Ramp_r(i) = Width(i)*delX*sum(q_Ramp_r);
    Q_Ramp_w(i) = Width(i)*delX*sum(q_Ramp_w);    
end

% Rho-Mu High-Speed method
x = [x1 x2 x3 x4 x5]; 
qr = [q_Tip_r q_Inlet_r q_Isolator_r q_Combustor_r q_Ramp_r]; 
Qr = [Q_Tip_r Q_Inlet_r+Q_Tip_r(end) Q_Isolator_r+Q_Inlet_r(end)+Q_Tip_r(end)  Q_Combustor_r+Q_Isolator_r(end)+Q_Inlet_r(end)+Q_Tip_r(end) Q_Ramp_r+Q_Combustor_r(end)+Q_Isolator_r(end)+Q_Inlet_r(end)+Q_Tip_r(end)];
qw = [q_Tip_w q_Inlet_w q_Isolator_w q_Combustor_w q_Ramp_w]; 
Qw = [Q_Tip_w Q_Inlet_w+Q_Tip_w(end) Q_Isolator_w+Q_Inlet_w(end)+Q_Tip_w(end)  Q_Combustor_w+Q_Isolator_w(end)+Q_Inlet_w(end)+Q_Tip_w(end) Q_Ramp_w+Q_Combustor_w(end)+Q_Isolator_w(end)+Q_Inlet_w(end)+Q_Tip_w(end)];
end
% figure 
% hold on
% plot(x,Qr-Qw,x,Qr,x,Qw)
% figure
% plot(x,qr-qw,x,qr,x,qw)