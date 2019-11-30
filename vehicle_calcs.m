%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   ERJ - Flight Calcs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Goals:
%   1. Need to figure out how far the vehicle can fly assuming
%      certain Thurst weight ratios
% 
%   2. Change out Isp values for different fuels

clear all;
close all;
clc;

%Reference

% 1. GHV - Liston Slides

q = [1000, 1500, 2000]*47.7; %Pa
%atmosphere fn takes ft., need to update every iteration
[T,P,Rho,dummy] = atmosphere4(Altitude,1);

T = T'*(5/9);   %K
P = P'.*47.77;  %Pa
Rho = Rho'*515.379; %kg/m3
R = 287;
gamma = 1.4;
Pr = 0.71;
r = Pr^(1/3);

% [1] Obtain Flight Parameters (Constant-q Trajectory)

for i = 1:length(q)
    a(:,i) = sqrt(gamma*R*T);
    U(:,i) = (sqrt(2*q(i)./Rho));
    M(:,i) = U(:,i)./a(:,i);
    P0(:,i) = P.*(1+(gamma-1)/2.*M(:,i).^2).^(gamma/(gamma-1));
    T0(:,i) = T.*(1+(gamma-1)/2.*M(:,i).^2);
    T0R(:,i) = 9/5.*T0(:,i);
    T0F(:,i) = T0R(:,i)-459.67;
    Tad(:,i) = T+r.*(T0(:,i)-T);
    TadR(:,i) = Tad(:,i).*9/5;
    TadF(:,i) = TadR(:,i) - 459.67;
end
 
% From GHV Vehicle:

% Mission Selection - Straight Run Out Accel to Upper Atmosphere Glide
% Specify Range to Accomplish Mission: Constantly accelerate upwards until
% run out of fuel

%Import the Aero Performance Data
Cd_Import

% Add Engine Performance Data (L&D) to Aero Data



%Run Mission Profiles for different Phi, Alpha, q, combinations & Graph??
