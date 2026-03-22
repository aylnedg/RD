function f = SpaceRendezvousAndDockingExampleObjectiveFcn( ...
    maneuverParameters, varargin)
%SPACERENDEZVOUSANDDOCKINGEXAMPLEOBJECTIVEFCN   Objective function
%   representing square of total delta-v magnitude that must be minimized
%   f = SpaceRendezvousAndDockingExampleObjectiveFcn(maneuverParameters, 
%      varargin) returns the objective function value, which is the square
%      of the sum of the magnitude of the two delta-vs.
%
%   Inputs:
%
%   maneuverparameters :8-by-1 vector of maneuver parameters. The first
%                       element represents the scaled time of first
%                       maneuver. Elements 2-4 represent the scaled delta-v
%                       for the first maneuver in the International
%                       Celestial Reference Frame (ICRF). Element 5
%                       represents the scaled time for the second maneuver
%                       counted from the first maneuver time. Elements 6-8
%                       represent the scaled delta-v for the second
%                       maneuver in ICRF.
%
%   varargin           :varargin represents the unused arguments supplied
%                       by fmincon and fsolve. These functions call
%                       objectiveFcn with 11 arguments, of which only the
%                       first argument, maneuverParameters, is used.

% Copyright 2023 The MathWorks, Inc.

% Extract the delta-v vectors.
deltaV1 = maneuverParameters(2:4);
deltaV2 = maneuverParameters(6:8);

% Calculate the objective function value.
f = (norm(deltaV1) + norm(deltaV2))^2;
end

