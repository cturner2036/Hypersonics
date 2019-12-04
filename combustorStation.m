function [To4,T4,Po4,P4,M4,pDrop] = combustorStation(To3, Po3, P3, T4, M3) 
maxIter = 5000; %number of segments
Pa = 2850; %Ambient Pressure (Pa)
To4 = 2400; %max stagnation temperature in scramjet (K)
Qr = 18875*1.05506/0.453592*1000; % J/kg heat of reaction of JP7 10% cracked, 50% gas
% To3 = 1646.7; %Isolator Exit Total Temperature (K)
% T3 = 1154.5; %Isolator Exit Static Temperature (K)
A3 = 0.0163135158; % Isolator Exit Area (m^2)
% Po3 = 1604141; %Isolator Exit Total Pressure (Pa)
% P3 = 462893.32; %Isolator Exit Static Pressure (Pa)
% M3 = 1.46; %Isolator Exit Mach
eta = 0.98; %scramjet combustor efficiency
ma = 22.661475; %mass flow of air (kg/s)
R = 287.05; %J/kg/K Gas Constant
% lc = 1; %length of combustor
[machi, Ti, Pi, rhoi, areai] = flowisentropic(1.4, 1.45, 'mach');
[machb, Tb, Pb, rhob, areab] = flowisentropic(1.28, 1.45, 'mach');
Cpb = CpPoly(To4.*Tb); %Combustor Cp (assuming minimal impact from fuel)
Cpi = CpPoly(To3.*Ti);
gammab = 1/(1 - R./Cpb);
gammai = 1/(1 - R./Cpi);
f = (Cpb.*To4./(Cpi.*To3) - 1)./( eta.*Qr./(Cpi.*To3) - Cpb.*To4./(Cpi.*To3));
mf = ma*f; %mass flow of fuel (kg/s)
mb = mf + ma; %total mass flow
% x = 0.01:0.01:1;
% theta = -5:0.1:10;
gamma = (gammai+gammab)/2; %average gamma of combustor
kc_arr = 0.5:0.0001:1.5; %combustor divergence parameter
for ii = 1:length(kc_arr)
    kc = kc_arr(ii);
    a = gamma*(kc-1)-kc;
    b = 2*kc - 1;
    alpha = 1 - 1/b;
    beta = 1/b;
    tau = 1;
    Mi = M3;
    M41 = 10;
    M42 = 1.3;
    Fo = To4/To3;
    iter = 1;
    while (abs(M41 - M42)>0.0001||iter==1)&&iter<maxIter
        M41 = M42;
        Fm = (a*M3^2 + b)/(a*M41^2 + b);
        Gm = M3^2/M41^2;
        Hm = (2 + (gamma - 1)*M41^2)/(2 + (gamma - 1)*M3^2);
        F1 = Fm^alpha*Gm^beta*Hm^tau;
        F = Fo - F1;
        dF1_dM4 = F1*( a*alpha/(a*M41^2 + b) + beta/M41^2 + (tau*(gamma - 1))/(2 + (gamma - 1)*M41^2));
        M42 = sqrt( M41^2 - F/dF1_dM4);
        iter = iter + 1;
    end
    ExitCondition(ii).Mach = M42;
    ExitCondition(ii).Ae = A3*(Fm^(-kc*(1/a + 1/b))*Gm^kc/b);
    ExitCondition(ii).P4 = P3*(Fm^((gamma*(kc - 1))/a));
    ExitCondition(ii).Po4 = Po3*(ExitCondition(ii).P4/P3)*Hm^(gamma/(gamma - 1));
    ExitCondition(ii).T4 = T3*(Fm^((b-1)/b)*Gm^(1/b));
    ExitCondition(ii).To4 = To3*(ExitCondition(ii).T4/T3)*Hm;
    if isreal(mb*sqrt(2*Cpb*ExitCondition(ii).To4*(1 - (Pa/ExitCondition(ii).Po4)^((gammab-1)/gammab))))&&iter<maxIter&&isreal(M42)
        ExitCondition(ii).Thrust = mb*sqrt(2*Cpb*ExitCondition(ii).To4*(1 - (Pa/ExitCondition(ii).Po4)^((gammab-1)/gammab)));
    else
        ExitCondition(ii).Thrust = NaN;
    end
end
% plot(kc_arr,[ExitCondition(:).Thrust]);
% xlabel('Divergence Factor kc')
% ylabel('Thrust (N)')
% grid on
[~,ind] = max([ExitCondition(:).Thrust]);
ExitCondition(ind);
To4 = ExitCondition(ind).To4;
T4 = ExitCondition(ind).T4;
Po4 = ExitCondition(ind).Po4;
P4 = ExitCondition(ind).P4;
M4 = ExitCondition(ind).Mach;
pDrop = abs(ExitCondition(ind).Po4 - Po3)/Po3;
​
function Cp = CpPoly(T)
    T = T.*1.8; %Convert to Rankine
    Cp = (3.2240709 + 2.5453604e-4.*T + 2.5222882e-7.*T.^2 -1.2013920e-10.*T.^3 + 1.7003460e-14.*T.^4).*1716.5;
    Cp = Cp.*0.1672;%J/kg/K Conversion
end
​
​
function F = gammaSolve(Po,To,A,mdot,M,gamma)
    F = A*Po/sqrt(To)*sqrt(gamma/287)*M*(1+(gamma-1)/2*M^2)^(-(gamma+1)/(2*(gamma-1))) - mdot;
end
​
​
​
​
end
