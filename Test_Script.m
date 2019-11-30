%Test Script

%Calculate Total Thrust from Engine
%   Ideal Nozzle
%[throat_height, throat_angle, cowl_height, body_width, step_size] = Plug_Nozzle(throat_angle, throat_height, cowl_height, body_width, step_size);
[throat_height,throat_angle,cowl_height,body_width,step_size] = Plug_Nozzle(67,0.2,1,1,100);
%   Truncated Nozzle
%[x,y] = Plug_Nozzle_Style2(AR,eta_b,throat_height);
[x,y,local_turn] = Plug_Nozzle_Style2(3.3839,0.05,throat_height,step_size);
%   Run Flow_Properties to calculate Thrust Values
%[Total_Thrust, Total_Lift] = Flow_Properties(step_size,local_turn,P_amb,T_exit,Pt_exit,throat_angle,throat_height,body_width,M_throat,y,x,alpha,DynamicPressure,m_dot,gamma);
[Total_Thrust, Total_Lift] = Flow_Properties(step_size,local_turn,1300,2400,1101325,throat_angle,throat_height,body_width,1.15,y,x,3,71820,22.675,1.4);

% Calculate total Thrust


%Calculate total Lift


%Update Vehicle Weight and Update Flight Values