function attitude = SpaceRendezvousAndDockingExampleCalculateAttitude( ...
    primaryAlignmentVector, primaryConstraintVector, ...
    secondaryAlignmentVector, secondaryConstraintVector)
%SpaceRendezvousAndDockingExampleCalculateAttitude   Calculate attitude
%   with respect to the International Celestial Reference Frame (ICRF)
%   based on the primary and secondary alignment and constraint vectors
%   attitude = SpaceRendezvousAndDockingExampleCalculateAttitude(
%      primaryAlignmentVector, primaryConstraintVector,
%      secondaryAlignmentVector, secondaryConstraintVector) calculates
%      returns scalar-first quaternion attitude with respect to ICRF based
%      on the specified primary and secondary alignment and constraint
%      vectors.
% 
%   Inputs:
%
%   primaryAlignmentVector    :3-by-1 vector defined in body frame that
%                              must be aligned with the primary constraint
%                              vector.
%   primaryConstraintVector   :3-by-1 vector defined in ICRF that the
%                              primary alignment vector must align with.
%   secondaryAlignmentVector  :3-by-1 vector defined in the body frame that
%                              must be aligned with the secondary
%                              constraint vector. Exact alignment is
%                              usually not possible for this vector.
%                              Therefore, the alignment error between the
%                              secondary alignment and constraint vectors
%                              must be minimized.
%   secondaryConstraintVector :3-by-1 vector defined in ICRF that the
%                              secondary constraint vector must aign with.
%
%   Outputs:
%
%   attitude :4-by-1 scalar-first quaternion representing the attitude with
%             respect to ICRF.

% Copyright 2023 The MathWorks, Inc.

% Convert the alignment and constraint vectors to unit vectors
primaryAlignmentVector = ...
    primaryAlignmentVector/norm(primaryAlignmentVector);
primaryConstraintVector = ...
    primaryConstraintVector/norm(primaryConstraintVector);
secondaryAlignmentVector = ...
    secondaryAlignmentVector/norm(secondaryAlignmentVector);
secondaryConstraintVector = ...
    secondaryConstraintVector/norm(secondaryConstraintVector);

% The desired alignment is achieved through two rotations. The first
% rotation aligns the primary alignment vector with the primary constraint
% vector. The desired rotation axis and angle of rotation can be calculated
% as follows.
if dot(primaryAlignmentVector, primaryConstraintVector) == 1
    % If primary constraint and alignment vectors are parallel, the first
    % rotation axis is taken as the primary constraint vector, and the
    % first rotation angle is 0.
    firstRotationAxis = primaryConstraintVector;
    firstRotationAngle = 0;
elseif dot(primaryAlignmentVector, primaryConstraintVector) == -1
    % If primary constraint and alignment vectors are anti-parallel, choose
    % a unit vector that is perpendicular to both. There are infinite such
    % vectors, just choose one and treat is as the first rotation axis. The
    % first rotation angle is pi.
    
    % Let the components of the primaryAlignmentVector vector are a1, a2, 
    % a3.
    a1 = primaryAlignmentVector(1);
    a2 = primaryAlignmentVector(2);
    a3 = primaryAlignmentVector(3);
    
    % Let the components of the first rotation axis be b1, b2, b3.
    if a1 ~= 0
        % If a1 is not equal to 0, choose b2 = 1, b3 = 1 if a2 ~ a3.
        % Otherwise, choose b2 = 1, b3 = -1.
        if a2 ~= a3
            b2 = 1;
            b3 = 1;
        else
            b2 = 1;
            b3 = -1;
        end
        % Solve for b3 such that the dot product between the primary
        % alignment vector and the first rotation axis is 0.
        b1 = (-(a2*b2) - (a3*b3))/a1;
    elseif a2 ~= 0
        % If a1 = 0, test a2. If it is not equal to 0, choose b1 = 0 and
        % b3 = 1. Solve for b2 such that the dot product between the 
        % primary alignment vector and the first rotation axis is 0.
        b1 = 0;
        b3 = 1;
        b2 = - (a3*b3/a2);
    else
        % If both a1 and a2 are zero, then the only requirement for the
        % first rotation axis to be perpendicular to the primary alignment
        % vector is that b3 = 0. Therefore, arbitrarily set b1 = 1 and b2 =
        % 0.
        b1 = 1;
        b2 = 0;
        b3 = 0;
    end
    % Construct the unit first rotation axis
    firstRotationAxis = [b1;b2;b3]/norm([b1;b2;b3]);
    
    % Set the first rotation angle to pi
    firstRotationAngle = pi;
else
    % If primary constraint and alignment vectors are neither parallel nor
    % anti-parallel, their cross product will result in a vector that is
    % perpendicular to both. The unit vector in the direction of the
    % product is the first rotation axis. The arc cosine of the dot product
    % of the primary alignment and constraint vectors will yield the first
    % rotation angle.
    firstRotationAxis = ...
        cross(primaryAlignmentVector, primaryConstraintVector);
    firstRotationAxis = firstRotationAxis/norm(firstRotationAxis);
    firstRotationAngle = acos(max(min( ...
        dot(primaryAlignmentVector, primaryConstraintVector), 1), -1));
end

% The second rotation attempts to minimize the angle between the secondary
% alignment vector and the secondary constraint vector.

% The second rotation is about the primary constraint vector.
% This ensures that the primary alignment vector after the first rotation 
% remains aligned with the primary constraint vector.
secondRotationAxis = primaryConstraintVector;

% Calculate the second rotation angle
if (dot(primaryConstraintVector, secondaryConstraintVector) == 1) || ...
        (dot(primaryConstraintVector, secondaryConstraintVector) == -1)
    % If the primary and secondary constraint vectors are parallel or
    % anti-parallel, the  the angle between the secondary alignment vector 
    % and secondary constraint vector will always remain at 90 degrees, 
    % regardless of the rotation angle. Therefore, the second rotation
    % angle is set to 0.
    secondRotationAngle = 0;
else
    % The rotation angle is calculated such that it minimizes the angle 
    % between the secondary alignment and constraint vectors.
    secondaryAlignmentVectorAfterFirstRotation = rotate( ...
        secondaryAlignmentVector, firstRotationAxis, ...
        firstRotationAngle);
    secondRotationAngle = rotationForAlignment( ...
        secondaryAlignmentVectorAfterFirstRotation, ...
        secondaryConstraintVector, secondRotationAxis);
end

% Calculated the body-fixed axes after the two rotations
xBody = [1;0;0];
xBody = rotate(xBody, firstRotationAxis, firstRotationAngle);
xBody = rotate(xBody, secondRotationAxis, secondRotationAngle);

yBody = [0;1;0];
yBody = rotate(yBody, firstRotationAxis, firstRotationAngle);
yBody = rotate(yBody, secondRotationAxis, secondRotationAngle);

zBody = [0;0;1];
zBody = rotate(zBody, firstRotationAxis, firstRotationAngle);
zBody = rotate(zBody, secondRotationAxis, secondRotationAngle);

% NED to Body transform
attitude = dcm2quat([xBody';yBody';zBody']);
end

function rotatedVector = rotate(vector, rotationAxis, rotationAngle)
%rotate   Helper function to rotate vector about rotationAxis by
%   rotationAngle
rotatedVector = vector*cos(rotationAngle) - ...
    cross(vector, rotationAxis*sin(rotationAngle)) + ...
    dot(vector, rotationAxis)*rotationAxis*(1 - cos(rotationAngle));
end

function theta = rotationForAlignment(vector1, vector2, rotationAxis)
%rotationForAlignment   Helper function to calculate the rotation angle
%   about rotationAxis that minimizes the angle between vector1 and vector2

% Extract the components of the rotation axis
lambda1 = rotationAxis(1);
lambda2 = rotationAxis(2);
lambda3 = rotationAxis(3);

% Extract the components of vector 1 (alignment vector)
a1 = vector1(1);
a2 = vector1(2);
a3 = vector1(3);

% Extract the components of vector 2 (constraint vector)
c1 = vector2(1);
c2 = vector2(2);
c3 = vector2(3);

% Analytically calculate the angle of rotation (theta) that minimizes the
% angle between vector1 and vector2. The derivation of the analytical
% expression is as follows:
%
% rotatedVector1 = rotatedVector = vector1*cos(theta) - ...
%     cross(vector, rotationAxis*sin(rotationAngle)) + ...
%     dot(vector, rotationAxis)*rotationAxis*(1 - cos(rotationAngle));
% projection = rotatedVector1.vector2
% 
% The projection needs to be maximized for alignment.
%
% d(projection)/drotationAngle = 0
%
% Solving the above equation for rotationAngle yields two solutions for 
% rotationAngle. The one that maximizes the projection is the desired 
% rotationAngle. Geometrically, these solutions correspond to two vectors
% that are parallel, but opposite to each other. Therefore, the projection
% will be positive for one, and negative for the other. Positive projection
% minimizes the angle between rotatedVector1 and vector2, while the
% negative one maximizes it.  Therefore, the rotationAngle that results in
% positive projection is chosen.

% The following expression is one of the analytical solutions of
% rotationAngle
theta = real(-log(-abs(a1*c1*lambda1^2 - a2*c2 - a3*c3 + ...
    a1*c2*lambda3*1i - a1*c3*lambda2*1i - a2*c1*lambda3*1i + ...
    a2*c3*lambda1*1i + a3*c1*lambda2*1i - a3*c2*lambda1*1i - a1*c1 + ...
    a2*c2*lambda2^2 + a3*c3*lambda3^2 + a1*c2*lambda1*lambda2 + ...
    a2*c1*lambda1*lambda2 + a1*c3*lambda1*lambda3 + ...
    a3*c1*lambda1*lambda3 + a2*c3*lambda2*lambda3 + ...
    a3*c2*lambda2*lambda3)/(a1*c1*lambda1^2 - a2*c2 - a3*c3 + ...
    a1*c2*lambda3*1i - a1*c3*lambda2*1i - a2*c1*lambda3*1i + ...
    a2*c3*lambda1*1i + a3*c1*lambda2*1i - a3*c2*lambda1*1i - a1*c1 + ...
    a2*c2*lambda2^2 + a3*c3*lambda3^2 + a1*c2*lambda1*lambda2 + ...
    a2*c1*lambda1*lambda2 + a1*c3*lambda1*lambda3 + ...
    a3*c1*lambda1*lambda3 + a2*c3*lambda2*lambda3 + ...
    a3*c2*lambda2*lambda3))*1i);

% Analytical expression for projection of vector rotatedVector1 on vector2
proj = c1*(a1*cos(theta) - a2*lambda3*sin(theta) + ...
    a3*lambda2*sin(theta) - lambda1*(cos(theta) - 1)*(a1*lambda1 + ...
    a2*lambda2 + a3*lambda3)) + c2*(a2*cos(theta) + ...
    a1*lambda3*sin(theta) - a3*lambda1*sin(theta) - lambda2*(cos(theta) ...
    - 1)*(a1*lambda1 + a2*lambda2 + a3*lambda3)) + c3*(a3*cos(theta) - ...
    a1*lambda2*sin(theta) + a2*lambda1*sin(theta) - lambda3*(cos(theta) ...
    - 1)*(a1*lambda1 + a2*lambda2 + a3*lambda3));

% If the projection is negative, it simply means that the solution used for
% rotation angle maximizes the angle between vector1 and vector2.
% Therefore, the other analytical solution will yield the desired rotation
% angle.
if proj < 0
    theta = real(-log(abs(a1*c1*lambda1^2 - a2*c2 - a3*c3 + ...
        a1*c2*lambda3*1i - a1*c3*lambda2*1i - a2*c1*lambda3*1i + ...
        a2*c3*lambda1*1i + a3*c1*lambda2*1i - a3*c2*lambda1*1i - ...
        a1*c1 + a2*c2*lambda2^2 + a3*c3*lambda3^2 + ...
        a1*c2*lambda1*lambda2 + a2*c1*lambda1*lambda2 + ...
        a1*c3*lambda1*lambda3 + a3*c1*lambda1*lambda3 + ...
        a2*c3*lambda2*lambda3 + ...
        a3*c2*lambda2*lambda3)/(a1*c1*lambda1^2 - a2*c2 - a3*c3 + ...
        a1*c2*lambda3*1i - a1*c3*lambda2*1i - a2*c1*lambda3*1i + ...
        a2*c3*lambda1*1i + a3*c1*lambda2*1i - a3*c2*lambda1*1i - ...
        a1*c1 + a2*c2*lambda2^2 + a3*c3*lambda3^2 + ...
        a1*c2*lambda1*lambda2 + a2*c1*lambda1*lambda2 + ...
        a1*c3*lambda1*lambda3 + a3*c1*lambda1*lambda3 + ...
        a2*c3*lambda2*lambda3 + a3*c2*lambda2*lambda3))*1i);
end
end
