function [throat_height, throat_angle, cowl_height, body_width, step_size, AR] = Plug_Nozzle(throat_angle, throat_height, cowl_height, body_width, step_size)
% Nozzle Contour, Approximate Method from Characteristic Line
%clc, clear all, close all
% Nozzle Contour, Approximate Method from Characteristic Line

% Read Me %
% This script calculates a perfectly expanded exhaust stream based
% on the thruster angle and and throat height design choices. This 
% will feed into the Rapid Method script in order to truncate the 
% nozzle at a specfic length.

% Input Dimensions

%throat_angle = 67;      % Throat Inclination
%throat_height = 0.2;    % Throat Height
%cowl_height = 1;        % Cowl Height to Centerline

% Secondary
%body_width = 1;         % Body Width
%step_size = 100;        % Model Fidelity

% throat_angle = 67;      % Throat Inclination
% throat_height = ;    % Throat Height
% cowl_height = 1;        % Cowl Height to Centerline
% 
% % Secondary
% body_width = 1;         % Body Width
% step_size = 100;        % Model Fidelity
gamma = 1.3;            % Ratio of Specific Heats

% Set  Matrices
x = zeros(1,step_size);
y = zeros(1,step_size);
mach_numbers = zeros(1,step_size);
mach_angles = zeros(1,step_size);
local_turn = zeros(1,step_size);
area_streamtube = zeros(1,step_size);
length_to_lip = zeros(1,step_size);

% Geometric Parameters
throat_area = throat_height*body_width;
turn_steps = linspace(0,throat_angle,step_size);

% Nozzle Contour Loop
for i = 1:step_size
    
    % Calculate Mach at all stations
    dummy_angle = 0;
    mach_numbers(i) = 1;
    while (dummy_angle < turn_steps(i))
        pm1 = sqrt((gamma+1)/(gamma-1));
        pm2 = (gamma-1)/(gamma+1);
        pm3 = mach_numbers(i)^2-1;           
        dummy_angle = pm1*atand(sqrt(pm2*pm3))-atand(sqrt(pm3));
        mach_numbers(i) = mach_numbers(i) + 0.00001;
    end
    
    % Calcualte Mach and Local Angles
    mach_angles(i) = asind(1/mach_numbers(i));
    local_turn(i) = throat_angle - turn_steps(i);
    
    % Calculate Streamtube Areas
    aa = (2+(gamma-1)*(mach_numbers(i)^2));
    bb = (2+(gamma-1)*(mach_numbers(1)^2));
    cc = (gamma+1)/(2*(gamma-1));
    area_streamtube(i) = throat_area * ( mach_numbers(1) / mach_numbers(i) ) * ( (aa / bb)^(cc) );
    
    % Calculate Contour
    length_to_lip(i) = area_streamtube(i)/sind(mach_angles(i));
    step_angle = mach_angles(i)+local_turn(i);
    x(i) = length_to_lip(i)*cosd(step_angle);
    y(i) = cowl_height - length_to_lip(i)*sind(step_angle);
   
end

x_hold = x(1);
y_hold = y(1);
for i = 1:step_size
    x(i) = x(i) - x_hold;
    y(i) = -1*(y(i) - y_hold);
end

% Find Expansion Ratio, Mexit, Check Cowl Height
AR = area_streamtube(end)/area_streamtube(1);
M_exit = mach_numbers(end);
y_height = y(1) - y(end);
x_length = x(end) - x(1);
d_cowl = sqrt( (x(1)-0)^2 + (y(1)-cowl_height)^2 );

% cout = ['Perfect Expansion Ramp Study. Expansion Ratio is: ',num2str(AR),...
%     '. Exit Mach is: ' ,num2str(M_exit),'. Only Ramp Height is: ', num2str(y_height),...
%     '. Total Length is: ',num2str(x_length),'.'];

%disp(cout)

% figure (1)
% plot(x,y)
% hold on
% plot(-x_hold,-1*(cowl_height-y_hold),'o')
% xlabel('Nozzle Length')
% ylabel('Nozzle Height')
% grid on
% axis equal
