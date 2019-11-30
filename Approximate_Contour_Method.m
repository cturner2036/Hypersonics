% Nozzle Contour, Approximate Method from Characteristic Line
clc, clear all, close all

% Input Dimensions
throat_angle = 25;      % Throat Inclination
throat_height = 0.2;    % Throat Height
cowl_height = 2;        % Cowl Height to Centerline
body_width = 1;         % Body Width
step_size = 100;        % Model Fidelity
gam = 1.4;              % Ratio of Specific Heats

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
        pm1 = sqrt((gam+1)/(gam-1));
        pm2 = (gam-1)/(gam+1);
        pm3 = mach_numbers(i)^2-1;           
        dummy_angle = pm1*atand(sqrt(pm2*pm3))-atand(sqrt(pm3));
        mach_numbers(i) = mach_numbers(i) + 0.00001;
    end
    
    % Calcualte Mach and Local Angles
    mach_angles(i) = asind(1/mach_numbers(i));
    local_turn(i) = throat_angle - turn_steps(i);
    
    % Calculate Streamtube Areas
    aa = (2+(gam-1)*(mach_numbers(i)^2));
    bb = (2+(gam-1)*(mach_numbers(1)^2));
    cc = (gam+1)/(2*(gam-1));
    area_streamtube(i) = throat_area * ( mach_numbers(1) / mach_numbers(i) ) * ( (aa / bb)^(cc) );
    
    % Calculate Contour
    length_to_lip(i) = area_streamtube(i)/sind(mach_angles(i));
    step_angle = mach_angles(i)+local_turn(i);
    x(i) = length_to_lip(i)*cosd(step_angle);
    y(i) = cowl_height - length_to_lip(i)*sind(step_angle);
    
end

figure (1)
plot(x,y)
xlabel('Nozzle Length')
ylabel('Nozzle Height')
grid on
