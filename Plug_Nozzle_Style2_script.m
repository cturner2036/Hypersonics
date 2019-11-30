% Nozzle Contour, Rapid Method for Plug Nozzle Design
clc, clear all

% Read me %
% This Rapid Method script creates the contour for a nozzle that is 
% truncated by a certain percentage. It requires the expansion ratio
% obtained from the Approximate Method and the throat height. This will
% feed into the flow properties script to calculate thrust and lift.

% Set Inputs
AR = 3.3839;                     % Expasion Area Ratio (Effective Angle)
eta_b = 0.1;                    % Truncation Parameter (0=PerfectExpansion)
throat_height = 0.2;             % Throat Height

% Secondary
step_size = 100;                 % Model Fidelity
gamma = 1.4;                     % Ratio Specific Heats
exit_height = AR*throat_height;  % Exit Height

% Set Matrices
turn_angle = zeros(1,step_size);
A = zeros(1,step_size);
l_nd = zeros(1,step_size);
r_nd = zeros(1,step_size);
x_nd = zeros(1,step_size);
y_nd = zeros(1,step_size);
l = zeros(1,step_size);
x = zeros(1,step_size);
y = zeros(1,step_size);
mach_numbers = zeros(1,step_size);
mach_angles = zeros(1,step_size);

% Exit Mach Solver (Copied from MST Thesis!)
f1 = 2/(gamma+1);
f2 = (gamma-1)/2;
f3 = (gamma+1)/(gamma-1);
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
pm1 = sqrt((gamma+1)/(gamma-1));
pm2 = (gamma-1)/(gamma+1);
pm3 = M^2-1;           
max_turn_angle = pm1*atand(sqrt(pm2*pm3))-atand(sqrt(pm3));
turn_steps = linspace(0,max_turn_angle,step_size);

% Contour Loop
for i = 1:step_size
    
    % Calculate Mach and PM Angles
    dummy_angle = 0;
    if 1~=i
       mach_numbers(i) = mach_numbers(i-1);
    else
       mach_numbers(i) = 1;
    end
    while (dummy_angle < turn_steps(i))
        pm1 = sqrt((gamma+1)/(gamma-1));
        pm2 = (gamma-1)/(gamma+1);
        pm3 = mach_numbers(i)^2-1;           
        dummy_angle = pm1*atand(sqrt(pm2*pm3))-atand(sqrt(pm3));
        mach_numbers(i) = mach_numbers(i) + 0.00001;
    end
    
    % Calcualte Mach and Alpha Angles
    mach_angles(i) = asind(1/mach_numbers(i));
    turn_angle(i) = turn_steps(end) - turn_steps(i) + mach_angles(i);
    
    % Calculate Streamtube Areas
    aa = (2+(gamma-1)*(mach_numbers(i)^2));
    bb = (2+(gamma-1)*(mach_numbers(1)^2));
    cc = (gamma+1)/(2*(gamma-1));
    A(i) = ( mach_numbers(1) / mach_numbers(i) ) * ( (aa / bb)^(cc) );
    
    % Calculate Non-Dimenional Parameters
    l_nd(i) = (1-sqrt(1-(A(i)*(1-(eta_b^2))...
        *mach_numbers(i)*(sind(turn_angle(i))/AR))))/sind(turn_angle(i));
    r_nd(i) = 1-(l_nd(i)*sind(turn_angle(i)));
    x_nd(i) = l_nd(i)*cosd(turn_angle(i));
    y_nd(i) = l_nd(i)*sind(turn_angle(i));
    
    % Calculate Coordinates
    if i == 1
        x_hold = x_nd(i)*exit_height;
        y_hold = y_nd(i)*exit_height;
    end
    l(i) = l_nd(i)*exit_height;
    x(i) = (x_nd(i)-x_nd(1))*exit_height;
    y(i) = (y_nd(i)-y_nd(1))*exit_height;
    
end

% Find Expansion Ratio
AR_check = A(end)/A(1);
M_exit_check = mach_numbers(end);
y_height = real(y(end))-y(1);
x_length = real(x(end))-x(1);

cout = ['Truncated Expansion Ramp Study. Expansion Ratio is: ',num2str(AR_check),...
    '. Exit Mach is: ' ,num2str(M_exit_check),'. Total Height is: ', num2str(y_height),...
    '. Total Length is: ',num2str(x_length),'.'];
disp(cout)

% Calcuate Cowl Position w.r.t. Ramp Origin
slope = (y(2)-y(1))/(x(2)-x(1));
slope = -1/slope;
n = 1000;
xc = linspace(0,0.4,n);
yc = slope.*xc;
for i = 1:n
    d = sqrt( xc(i)^2 + yc(i)^2);
    if d > throat_height
        xnew = xc(i);
        ynew = yc(i);
        break;
    end
end

local_turn = flip(turn_steps);

figure(1)
plot(x,y)
hold on
plot(xnew,ynew,'o')
xlabel('Nozzle Length')
ylabel('Nozzle Height')
grid on
axis equal



