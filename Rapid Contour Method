% Nozzle Contour, Rapid Method for Plug Nozzle Design
clc, clear all, close all

% Set Inputs
AR = 2;           % Expasion Area Ratio
eta_b = 0.01;     % Truncation Parameter (0=PerfectExpansion)
t_diam = 0.2;     % Throat Height (UseWhateverUnits)
step_size = 100;  % Model Fidelity
gam = 1.4;        % Ratio Specific Heats
 
% Derived Geometry
A_t = pi*((t_diam/2)^2);    % Throat Area - Assumed Circle
A_e = AR*A_t;               % Exit Area
re = sqrt(A_e/pi);          % Exit Radius

% Set Matrices
alpha = zeros(1,step_size);
A = zeros(1,step_size);
l_nd = zeros(1,step_size);
r_nd = zeros(1,step_size);
x_nd = zeros(1,step_size);
y_nd = zeros(1,step_size);
lr = zeros(1,step_size);
r = zeros(1,step_size);
xr = zeros(1,step_size);
yr = zeros(1,step_size);
mach_numbers = zeros(1,step_size);
mach_angles = zeros(1,step_size);

% Exit Mach Solver (Copied from MST Thesis!)
f1 = 2/(gam+1);
f2 = (gam-1)/2;
f3 = (gam+1)/(gam-1);
M0 = 4;
M = M0;
Function = (1/(M^2))*((f1*(1+(f2*(M^2))))^f3)-(AR^2);
Tolerance = 0.0001;
maxits = 200;
J = abs(Function);
i = 1;
while J>Tolerance && i<=maxits
    Function = (1/(M^2))*((f1*(1+(f2*(M^2))))^(f3))-(AR^2);
    dFunction = (f3/(M^2))*((f1*(1+(f2*(M^2))))^(f3-1))*2*f2*f1*M + (-2/(M^3))*((f1*(1+f2*(M^2)))^(f3));
    y = M-(Function/dFunction);
    M = y;
    Function = (1/(M^2))*((f1*(1+(f2*(M^2))))^(f3))-(AR^2);
    J = abs(Function);
    i = i+1;
end

% Use P-M to find Max Turn Angle from Exit Mach
pm1 = sqrt((gam+1)/(gam-1));
pm2 = (gam-1)/(gam+1);
pm3 = M^2-1;           
max_turn_angle = pm1*atand(sqrt(pm2*pm3))-atand(sqrt(pm3));
turn_steps2 = linspace(0,max_turn_angle,step_size);

% Contour Loop
for i = 1:step_size
    
    % Calculate Mach and PM Angles
    dummy_angle = 0;
    if 1~=i
       mach_numbers(i) = mach_numbers(i-1);
    else
       mach_numbers(i) = 1;
    end
    while (dummy_angle < turn_steps2(i))
        pm1 = sqrt((gam+1)/(gam-1));
        pm2 = (gam-1)/(gam+1);
        pm3 = mach_numbers(i)^2-1;           
        dummy_angle = pm1*atand(sqrt(pm2*pm3))-atand(sqrt(pm3));
        mach_numbers(i) = mach_numbers(i) + 0.00001;
    end
    
    % Calcualte Mach and Alpha Angles
    mach_angles(i) = asind(1/mach_numbers(i));
    alpha(i) = turn_steps2(end) - turn_steps2(i) + mach_angles(i);
    
    % Calculate Streamtube Areas
    aa = (2+(gam-1)*(mach_numbers(i)^2));
    bb = (2+(gam-1)*(mach_numbers(1)^2));
    cc = (gam+1)/(2*(gam-1));
    A(i) = ( mach_numbers(1) / mach_numbers(i) ) * ( (aa / bb)^(cc) );
    
    % Calculate Non-Dimenional Parameters
    l_nd(i) = (1-sqrt(1-(A(i)*(1-(eta_b^2))...
        *mach_numbers(i)*(sind(alpha(i))/AR))))/sind(alpha(i));
    r_nd(i) = 1-(l_nd(i)*sind(alpha(i)));
    x_nd(i) = l_nd(i)*cosd(alpha(i));
    y_nd(i) = l_nd(i)*sind(alpha(i));
    
    % Calculate Coordinates
    lr(i) = l_nd(i)*re;
    r(i) = r_nd(i)*re;
    xr(i) = x_nd(i)*re*10;
    yr(i) = y_nd(i)*re*10;
    
end

figure(1)
plot(xr,yr)
xlabel('Nozzle Length')
ylabel('Nozzle Height')
grid on



