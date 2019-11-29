function [realtime_weight] = actualWeight(full_weight, mf_dot, timestep)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function is meant to be implemented in a timestep loop
%Variables:
% - full weight: mass of vehicle at flight initiation
% - mf_dot: rate of fuel leaving the vehicle
%%%
%Output:
% - current weight of aircraft system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
realtime_weight = full_weight-mf_dot*timestep;