%% HW1 homework Problem #2 SCRamjet Solution by Joshua E. McMillin
%% Free Stream (Station 0)
hgeom(:,1)=0:100:299200;
[T0,P0,rho,Hgeop]=atmosphere4(hgeom,1);
gamma=1.4;
q=1500;
CPR = 10;
ExitLimit = 4459.67; % R = 3500F
M0=((2./(gamma*P0))*q).^0.5;
figure
plot(M0,hgeom)
APlot=gca;
grid on;
xlim(APlot,[0 25]);
title(APlot,'Altitude vs. Mach for q=1500psf');
xlabel(APlot,'Mach');
ylabel(APlot,'Geometric Altitude (ft)');
tau_r =(1+0.5*(gamma-1)*(M0.*M0)); 
T0t = T0 .* tau_r;
pi_r = tau_r.^(gamma/(gamma-1));
P0t =P0 .* pi_r;

%% Station 2
% inlet recovery MIL-E-5008B
pi_i = ones(size(M0)); % M.le.1 pi_i = 1.0 
pi_i(M0>1.0) = 1.0-0.075.*(M0(M0>1)-1).^1.35; % 1.lt.M.ge.5 pi_i = 1.0-0.075*(M0-1)^1.35
pi_i(M0>5.0)= 800./(935 + M0(M0>5.0).^4); % M.gt.5 pi_i = 800/(935+ M0^4)

P2t = P0.*pi_r.*pi_i;
T2t = T0 .* tau_r;

%% Station 4 %Combustion Exit
P4t = P0.*pi_r.*pi_i; % = P2t
T4t = ExitLimit;
Cp = 0.3; %BTU/lbm-degF Specific heat of air
dh_c = 18500; %BTU/lbm Heat of combustion per unit mass of fuel
% dm_f * dh_c = dm_air * Cp(T4t-T3t)
fc = Cp.*(T4t-T2t) ./ dh_c; %Fuel to air Ratio

%% Station 9
P9t = P0 .* pi_r .* pi_i; % P9t=P4t
%P9=P0
T9t = T4t; %Exit Limit
T9 = T9t .*((P0./P9t).^((gamma-1)/gamma));
% P9t/P9 = (1+0.5*(gamma-1)*M9^2)^(gamma/(gamma-1))
M9 = ((2/(gamma-1))*((P9t./P0).^((gamma-1)/gamma)-1)).^0.5;

%% Isp
% F = m(1+fc)U9 - mU0
% F/(m*fc) = (U0/fc) * ((1+fc) U9/U0 - 1)
% Isp = M0/fc *(gamma * R *T0)^.5 * ((1+fc)(M9/M0)(T9/T0)^0.5 - 1)
%Isp = M0/fc *(gamma * P0/rho_0)^.5 * ((1+fc)(M9/M0)(T9/T0)^0.5 - 1)
Isp = ((1+fc).*(M9./M0).*(T9./T0).^0.5 - 1).* (M0./fc) .*((gamma.*P0./rho).^0.5)./32.174; % F/m_air
figure;
plot(M0(and(M0>3.5,M0<7)),Isp(and(M0>3.5,M0<7)))
IPlot=gca;
grid on;
%xlim(IPlot,[1.7 4]);
title(IPlot,'Isp vs. Mach for q=1500psf');
xlabel(IPlot,'Mach');
ylabel(IPlot,'{I}_{sp} ({seconds})');
SCRJet = [M0(and(M0>3.5,M0<7)) Isp(and(M0>3.5,M0<7))];
save('SCRJet.mat','SCRJet')
