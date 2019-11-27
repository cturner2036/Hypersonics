% Calculate Isentropic Ramp Flow Propertues
% Works with Approximate Method Values, WIP for Rapid Contour Script
% Run Approx Method Value First, then this script

%% Updates
%%% INCLUDED LIFT AND DRAG AS FUNCTION OF ALPHA
%%% Calculated volume from dimensions of the nozzle!!!
%%% Calculated center of pressure for nozzle centerbody

%% Inputs (Arbitrary)
% Exit Flow Properties
M_throat = 1.15;    % Throat Mach
mdot = 50/2.205;    % Mass Flow Rate (kg/s)
R = 287;            % Specific Gas Constant (J/kg-K)
gamma = 1.4;        % Ratio Specific Heats
Q = 1500*47.88;     % Dynamic Pressure (Pa)
alpha = 3;          % Angle of Attack

% Importanat, will change with altitude!!!!
P_amb = 1300;     % Ambient-Inf Pressure (Pa)

% Static and Stagnation, Temperature and Pressure (K and Pa)
T_exit = 2400;      
Tt_exit = T_exit*(1 + (gamma-1)/2*M_throat^2);
Pt_exit = 1101325;
P_exit = Pt_exit*(1 + (gamma-1)/2*M_throat^2)^(-gamma/(gamma-1));

% Also requires: throat_angle, throat_area, local_turn, step_size, body_width

%-----------------------------------------------------------%
%% Initialize Parameters
% Set Matrices
Mach_Numbers = zeros(1,step_size);
Static_Temperature = zeros(1,step_size);
Stagnation_Temperature = zeros(1,step_size);
Static_Pressure = zeros(1,step_size);
Stagnation_Pressure = zeros(1,step_size);
Cof_Pressure = zeros(1,step_size);
Cof_AxialForce = zeros(1,step_size);
Cof_NormalForce = zeros(1,step_size);
Cof_Lift = zeros(1,step_size);
Cof_Thrust = zeros(1,step_size);

% Initialize Combustor Exit Properties
Mach_Numbers(1) = M_throat;
Static_Temperature(1) = T_exit;
Stagnation_Temperature(1) = Tt_exit;
Static_Pressure(1) = P_exit;
Stagnation_Pressure(1) = Pt_exit;

% Geometry Values
Base_height = y(end);
Base_area = Base_height*body_width;

% Exit Combustor Velocity
a = sqrt(gamma*R*T_exit);
v_exit = M_throat*a;

%% Calcualtes Mach Loop
% Required if exit Mach != 1
for i = 1:step_size-1
    
    % P-M Relations for v(M1)
    pm1 = sqrt((gamma+1)/(gamma-1));
    pm2 = (gamma-1)/(gamma+1);
    pm3 = Mach_Numbers(i)^2-1;           
    v1 = pm1*atand(sqrt(pm2*pm3))-atand(sqrt(pm3));
    
    % Local Turn Angle
    % Theta = v(M2) - v(M1)
    Theta = local_turn(i) - local_turn(i+1);
    v2 = v1 + Theta;
    
    % P-M Relations for Mach after Isentropic Turn
    Mach_Numbers(i+1)=0;
    dummy_angle = 0;
    while (dummy_angle < v2)
        pm1 = sqrt((gamma+1)/(gamma-1));
        pm2 = (gamma-1)/(gamma+1);
        pm3 = Mach_Numbers(i+1)^2-1;           
        dummy_angle = pm1*atand(sqrt(pm2*pm3))-atand(sqrt(pm3));
        Mach_Numbers(i+1) = Mach_Numbers(i+1) + 0.00001;
        % Provides three decimal place accuracy comparable w/ onlince calculators
    end
    
end


%% Flow Properties Loop
% Assuming Isentropic Expansion, Adiabatic Wall
for i = 1:step_size-1
    
    % Isentropic Relations
    TopVal = 1 + (gamma-1)/2*Mach_Numbers(i)^2;
    BotVal = 1 + (gamma-1)/2*Mach_Numbers(i+1)^2;
    exp1 = gamma/(gamma-1);
    
    % Static Pressure
    Static_Pressure(i+1) = Static_Pressure(i)*(TopVal/BotVal)^(exp1);
    
    % Stagnation Pressure
    Stagnation_Pressure(i+1) = Static_Pressure(i+1)*TopVal^(exp1);
    
    % Static Temperature
    Static_Temperature(i+1) = Static_Temperature(i)*(TopVal/BotVal);
    
    % Stagnation Temperature
    Stagnation_Temperature(i+1) = Static_Temperature(i+1)*(TopVal);
    
end 

%% Centerbody Coefficient Loop
for i = 1:step_size
    
    % Pressure Coefficinent
    Cof_Pressure(i) = (Static_Pressure(i)-P_amb)/Q;

    % Axial Force Coefficient
    Cof_AxialForce(i) = Cof_Pressure(i)*sind(local_turn(i));
    
    % Normal Force Coefficient
    Cof_NormalForce(i) = Cof_Pressure(i)*cosd(local_turn(i));
    
    % Lift Coefficient
    Cof_Lift(i) = Cof_NormalForce(i)*cosd(alpha) + Cof_AxialForce(i)*sind(alpha);
    
    % Thrust Coefficient
    Cof_Thrust(i) = Cof_AxialForce(i)*cosd(alpha) - Cof_NormalForce(i)*sind(alpha);
    
end

%% Calculate Thrust and Lift Components
% Jet Thrust and Lift
T_jet = (mdot*v_exit+(Static_Pressure(1)-P_amb)*(throat_area))*cosd(throat_angle-alpha);
L_jet = (mdot*v_exit+(Static_Pressure(1)-P_amb)*(throat_area))*sind(throat_angle-alpha);
% Base Thrust and Lift
T_base = Base_area*(Static_Pressure(end)-P_amb)*cosd(alpha);
L_base = Base_area*(Static_Pressure(end)-P_amb)*sind(alpha);

L_pressure = 0;
L_pressure2 = 0;
T_pressure = 0;
T_pressure2  = 0;
% Centerbody Thrust
for i = 1:step_size-1
       
    % First Method
    % Axial Cross Sectional Slice
    dAy = body_width*(y(i)-y(i+1));
    % Total Pressure on Slice
    T_dP = dAy*Q*Cof_Pressure(i)*cosd(alpha);
    % Centerbody Thrust
    T_pressure = T_pressure + T_dP;
    % Calculate Lift w/ First Method
    dAx = body_width*abs( (abs(x(i))-abs(x(i+1))) );
    L_dP = dAx*Q*Cof_Pressure(i)*cosd(alpha);
    L_pressure = L_pressure + L_dP;
    
    % Second Method
    % Total Area of Section
    Area2 = body_width*(sqrt( (x(i)-x(i+1))^2 + (y(i)-y(i+1))^2 ));
    % Axial Pressure on Slice
    T_dP2 = Area2*Q*Cof_Thrust(i);
    L_dP2 = Area2*Q*Cof_Lift(i);
    % Centerbody Thrust
    T_pressure2 = T_pressure2 + T_dP2;
    L_pressure2 = L_pressure2 + L_dP2;

end

% Combined Thrust and Lift Term (kN)
Total_Thrust = (T_pressure + T_jet + T_base)/1000;
Total_Lift = (L_pressure + L_jet + L_base)/1000;


%% Center of Pressure - Centerbody
% Numerically Integrate S x*P dx
xPx = trapz(x,x.*Static_Pressure);
% Numerically Integrate S P dx
Px = trapz(x,Static_Pressure);
% Find Center of Pressure, Useful for Moment Calculation
CoPx = xPx./Px;

%% Calculate Volume of Nozzle
Area = trapz(x,y);
Volume = Area*body_width;

figure (2)
plot(x,Static_Pressure)
hold on
plot(x,Stagnation_Pressure)
hold off
xlabel('Axial Distance')
ylabel('Pressure (Pa)')
legend("Static Pressure","Stagnation Pressure")
grid on

figure (3)
plot(x,Static_Temperature)
hold on
plot(x,Stagnation_Temperature)
hold off
xlabel('Axial Distance')
ylabel('Temperature (K)')
legend("Static Temperature","Stagnation Temperature")
grid on

figure (4)
plot(x,Cof_NormalForce)
hold on
plot(x,Cof_AxialForce)
hold off
xlabel('Axial Distance')
ylabel('Force Coefficients')
legend("Normal Force Coefficient","Axial Force Coefficinet")
grid on


