# Hypersonics
Purdue AAE 537 - Hypersonics Final Project Code
%% Define constants

%% Flight Simulator - Need to settle on when we update values such as velocity terms, AoA and mass properties
  % calculate flight conditions
  X = x0 + dx*dt + 0.5*dx2 * dt^2;
  h = h0 + dh*dt + 0.5*dh2 * dt^2;
  dx = dx + dx2*dt;
  dh = dh + dh2*dt;
  v = sqrt(dx^2, dh^2)
  
  % Calculate ambient conditions
  [temp,press,rho,~]=atmosphere4(h,GeometricFlag)
  Q = 0.5*rho*v^2 % watch the units here rho i sin slug/ft^2
  M = v/sqrt(gamma*R*T); Flight Mach number
  
  %convert to SI units???
  
%% vehicle simulator
  % is angle of attack 0 or are we trying to controlsomething?
    % What about control surfaces do we add drag to balance out moments induced by lift forces?
  % update mass properties for fuel consumption
  
  If %fuel is above bingo level
    %calculate drag
      %from vehicle definition
    %Calculate thrust
      %get inlet conditions from GHV
      %Get combustor outlet conditions and fuel flow from combustor model
      % calculate nozzle thrust & lift
   end
   
   %sum axial force and drag = dx^2
   %sum lift forces and drag = dh^2
   
   %Report conditions and store time history
   % repeat
    


