function varargout = SpaceRendezvousAndDockingExampleConstraintFcn( ...
    maneuverParemeters, epoch, initialChaserOscullatingElements, ...
    initialTargetOscullatingElements, finalRelativeRadialDiscance, ...
    finalRelativeVelocityICRF, minimumChaserPeriapsisBetweenManeuvers, ...
    needToReturnInequalityConstraint, scalePos, scaleVel, scaleTime)
%SPACERENDEZVOUSANDDOCKINGEXAMPLECONSTRAINTFCN  Constraint function
%   [c,ceq] = SpaceRendezvousAndDockingExampleConstraintFcn(
%      maneuverParemeters, epoch, initialChaserOscullatingElements, 
%      initialTargetOscullatingElements, finalRelativeRadialDiscance, 
%      finalRelativeVelocityICRF, minimumChaserPeriapsisBetweenManeuvers,
%      needToReturnInequalityConstraint, scalePos, scaleVel, scaleTime)
%      returns the inequality and equality constrints (c and ceq) when
%      needToReturnInequalityConstraint is true. These values are
%      calculated by extracting and applying the maneuver parameters on the
%      chaser based on its initial oscullating elements, and calculating
%      the states of the target after the second maneuver on the chaser.
%
%   ceq = SpaceRendezvousAndDockingExampleConstraintFcn(maneuverParemeters,
%      epoch, initialChaserOscullatingElements,
%      initialTargetOscullatingElements, finalRelativeRadialDiscance, 
%      finalRelativeVelocityICRF, minimumChaserPeriapsisBetweenManeuvers,
%      needToReturnInequalityConstraint, scalePos, scaleVel, scaleTime)
%      returns the inequality and equality constrints (c and ceq) when
%      needToReturnInequalityConstraint is false.
%
%   Inputs:
%
%   maneuverparameters                     :8-by-1 vector of maneuver
%                                           parameters. The first element
%                                           represents the scaled time of
%                                           first maneuver. Elements 2-4
%                                           represent the scaled delta-v
%                                           for the first maneuver in the
%                                           International Celestial
%                                           Reference Frame (ICRF). Element
%                                           5 represents the scaled time
%                                           for the second maneuver counted
%                                           from the first maneuver time.
%                                           Elements 6-8 represent the
%                                           scaled delta-v for the second
%                                           maneuver in ICRF.
%   epoch                                  :Mission start time, specified
%                                           as a scalar datetime object.
%   initialChaserOscullatingElements       :Initial oscullating elements of
%                                           the chaser, specified as a
%                                           6-by-1 vector consisting of
%                                           semimajor axis (m),
%                                           eccentricity, inclination
%                                           (deg), right ascension of
%                                           ascending node (deg), argument
%                                           of periapsis (deg) and true
%                                           anomaly(deg).
%   initialTargetOscullatingElements       :Initial oscullating elements of
%                                           the target, specified as a
%                                           6-by-1 vector consisting of
%                                           semimajor axis (m),
%                                           eccentricity, inclination
%                                           (deg), right ascension of
%                                           ascending node (deg), argument
%                                           of periapsis (deg) and true
%                                           anomaly(deg).
%   finalRelativeRadialDiscance            :Relative distance of chaser
%                                           with respect to the target
%                                           after the second maneuver along
%                                           the radial position vector of
%                                           the target, specified as a
%                                           scalar in m.
%   finalRelativeVelocityICRF              :Final relative velocity of
%                                           chaser with respect to target
%                                           in ICRF, specified as a 3-by-1
%                                           vector in m/s.
%   minimumChaserPeriapsisBetweenManeuvers :Minimum periapsis of chaser
%                                           inbetween maneuvers, specified
%                                           as a scalar in m.
%   needToReturnInequalityConstraint       :Flag to determine whether
%                                           inequality constraint must be
%                                           returned, specified as a
%                                           logical scalar.
%   scalePos                               :Position scale, specified as a
%                                           scalar in m.
%   scaleVel                               :Velocity scale, specified as a
%                                           scalar in m/s.
%   scaleTime                              :Time scale, specified as a
%                                           scalar in s.
%
%   Outputs:
%
%   c   :Inequality constraint on chaser periapsis between maneuvers,
%        specified as a scalar. The value returned equals:
%        - (chaserPeriapsisBetweenManeuvers 
%        - minimumChaserPeriapsisBetweenManeuvers)/scalePos
%        c is returned only when needToReturnInequalityConstraint is true.
%   ceq :2-by-1 equality constraint vector. The first element represents
%        equality constraint on final ICRF position of chaser given by:
%        (finalChaserICRFPosition -
%        (finalTargetICRFPosition*finalRelativeRadialDiscance/
%        norm(finalTargetICRFPosition)))/scalePos
%        The second element represents equality constraint on final ICRF
%        velocity of chaser given by:
%        (finalChaserICRFVelocity - finalTargetICRFVelocity -
%        finalRelativeVelocityICRF)/scaleVel

% Copyright 2023 The MathWorks, Inc.

% Calculate the initial position and velocity of the chaser from the
% oscullating elements.
[positionChaser, velocityChaser] = keplerian2ijk( ...
    initialChaserOscullatingElements(1), ...
    initialChaserOscullatingElements(2), ...
    initialChaserOscullatingElements(3), ...
    initialChaserOscullatingElements(4), ...
    initialChaserOscullatingElements(5), ...
    initialChaserOscullatingElements(6));

% Calculate the burn time for the first maneuver.
tManeuver = maneuverParemeters(1)*scaleTime;
firstBurnTime = epoch + seconds(tManeuver);

% Extract the first delta-V vector.
deltaV = maneuverParemeters(2:4)*scaleVel;

% Calculate chaser states after first maneuver.
[positionChaser, velocityChaser] = propagateOrbit( ...
    firstBurnTime, ...
    positionChaser, ...
    velocityChaser, ...
    PropModel = 'two-body-keplerian', ...
    Epoch = epoch);
velocityChaser = velocityChaser + deltaV;

% Calculate chaser periapsis after the first maneuver.
[semiMajorAxis, eccentricity] = ...
    ijk2keplerian(positionChaser, velocityChaser);
chaserPeriapsisBetweenManeuvers = semiMajorAxis*(1 - eccentricity);

% Calculate the second burn time.
tManeuver = maneuverParemeters(5)*scaleTime;
secondBurnTime = firstBurnTime + seconds(tManeuver);

% Extract the second delta-V vector.
deltaV = maneuverParemeters(6:8)*scaleVel;

% Calculate chaser states after second maneuver.
[positionChaser, velocityChaser] = propagateOrbit( ...
    secondBurnTime, ...
    positionChaser, ...
    velocityChaser, ...
    PropModel = 'two-body-keplerian', ...
    Epoch = firstBurnTime);
velocityChaser = velocityChaser + deltaV;

% Calculate the states of the target spacecraft at final time.
[positionTarget, velocityTarget] = propagateOrbit( ...
    secondBurnTime, ...
    initialTargetOscullatingElements(1), ...
    initialTargetOscullatingElements(2), ...
    initialTargetOscullatingElements(3), ...
    initialTargetOscullatingElements(4), ...
    initialTargetOscullatingElements(5), ...
    initialTargetOscullatingElements(6), ...
    PropModel = 'two-body-keplerian', ...
    Epoch = epoch);

% Position error of chaser after the second maneuver.
desiredPosition = positionTarget + ...
    finalRelativeRadialDiscance*positionTarget/norm(positionTarget);
positionError = positionChaser - desiredPosition;

% Velocity error of chaser after the second maneuver.
velocityError = ...
    velocityChaser - velocityTarget - finalRelativeVelocityICRF;

% Equality constraint.
ceq = [positionError/scalePos; ...
    velocityError/scaleVel];

% Inequality constraint.
c = - (chaserPeriapsisBetweenManeuvers -  ...
    minimumChaserPeriapsisBetweenManeuvers)/scalePos;

% Prepare the outputs.
if needToReturnInequalityConstraint
    varargout = {c, ceq};
else
    varargout = {ceq};
end
end
