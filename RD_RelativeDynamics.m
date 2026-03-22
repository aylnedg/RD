function [xdot, aux] = RD_RelativeDynamics(t, x, u_force, u_torque, params)
%RD_RELATIVEDYNAMICS  Relative translational and 6-DOF rotational dynamics
%   for all R&D mission phases.
%
%   [xdot, aux] = RD_RelativeDynamics(t, x, u_force, u_torque, params)
%
%   Implements three dynamics models selectable via params.mode:
%
%   1) 'CW'    – Clohessy-Wiltshire linearised relative motion (Phases 3,4)
%      State x = [dx; dy; dz; dxdot; dydot; dzdot]  (6x1)
%      LVLH frame: x=radial(outward), y=along-track(velocity), z=cross-track
%
%   2) '6DOF'  – Full nonlinear 6-DOF chaser dynamics (Phases 5,6)
%      State x = [r_ECI(3); v_ECI(3); q_bi(4); omega_b(3)]  (13x1)
%
%   3) 'keplerian' – Two-body relative motion via absolute Keplerian propagation
%      (For initialisation / far-range checks)
%      State x = [r_ECI(3); v_ECI(3)]  (6x1)
%
%   Inputs:
%     t         scalar      simulation time [s]
%     x         [Nx1]       state vector (see above)
%     u_force   [3x1] N     external force in body / LVLH frame
%     u_torque  [3x1] Nm    external torque in body frame (6DOF only)
%     params    struct       model parameters (see below)
%
%   params fields:
%     .mode          'CW' | '6DOF' | 'keplerian'
%     .n             mean motion [rad/s]   (required for CW)
%     .mu            gravitational param  [m³/s²]  (required for keplerian/6DOF)
%     .mass_kg       chaser mass [kg]
%     .Inertia_kgm2  [3x3] inertia tensor (required for 6DOF)
%     .J2_enable     include J2 oblateness perturbation (bool)
%     .drag_enable   include atmospheric drag (bool)
%     .CD            drag coefficient (required if drag_enable)
%     .A_m2          cross-sectional area [m²] (drag)
%     .rho_kgm3      atmospheric density [kg/m³] (drag, approx. constant)
%
% Reference:
%   Fehse, W. (2003). Automated Rendezvous and Docking of Spacecraft.
%   Cambridge University Press. Ch. 3, 4.
%   Clohessy & Wiltshire (1960). Terminal guidance system for satellite
%   rendezvous. J. of the Aerospace Sciences, 27(9), 653-658.
%
% See also: RD_GuidanceControl, RD_NavigationFilter, RD_RunAllPhases

% Copyright 2024 - Extended R&D Simulation

% Default parameters
mu   = getparam(params, 'mu',      3.986004418e14); % m³/s²  Earth
n    = getparam(params, 'n',       1.1e-3);         % rad/s  (ISS orbit ~400 km)
mass = getparam(params, 'mass_kg', 8000);           % kg

xdot = zeros(size(x));
aux  = struct('acc_grav', zeros(3,1), 'acc_drag', zeros(3,1), ...
              'acc_J2', zeros(3,1));

mode = getparam(params, 'mode', 'CW');

switch mode

    % ===================================================================
    case 'CW'
    % -------------------------------------------------------------------
    % Clohessy-Wiltshire equations (linearised HCW equations)
    % Valid for: small relative distances (< ~5 km), near-circular target orbit
    %
    % State: x = [dx; dy; dz; dxdot; dydot; dzdot]
    %   dx    = radial      separation (x-axis outward)
    %   dy    = along-track separation (y-axis in velocity dir)
    %   dz    = cross-track separation (z-axis = h/|h|)
    %   dxdot, dydot, dzdot = time derivatives in LVLH frame
    %
    % CW equations (Fehse Eq. 3.38):
    %   ddx   = 3n²dx + 2n*dydot + fx/m
    %   ddy   = -2n*dxdot       + fy/m
    %   ddz   = -n²*dz           + fz/m
    % -------------------------------------------------------------------

        if numel(x) ~= 6
            error('CW mode requires 6-element state vector');
        end

        dx   = x(1); dy  = x(2); dz   = x(3);
        dxd  = x(4); dyd = x(5); dzd  = x(6);

        % External force per unit mass in LVLH (from RCS)
        fx = u_force(1) / mass;
        fy = u_force(2) / mass;
        fz = u_force(3) / mass;

        % CW dynamics
        xdot(1) = dxd;
        xdot(2) = dyd;
        xdot(3) = dzd;
        xdot(4) = 3*n^2*dx + 2*n*dyd + fx;   % radial
        xdot(5) = -2*n*dxd            + fy;   % along-track
        xdot(6) = -n^2*dz             + fz;   % cross-track

        aux.acc_grav = [3*n^2*dx; 0; -n^2*dz];

        % Optional: tidal / J2 perturbation correction term
        if getparam(params, 'J2_enable', false)
            R_E = 6.371e6; % m
            J2  = 1.08263e-3;
            % Simplified first-order J2 correction in LVLH (Schaub & Junkins)
            r0  = getparam(params, 'target_radius_m', 6.771e6);
            delta_J2 = (3/2)*J2*(mu*R_E^2)/(r0^5) * ...
                       [-4*dx; 2*dy; -dz];  % approx. in-plane dominant
            xdot(4:6) = xdot(4:6) + delta_J2/mass;
            aux.acc_J2 = delta_J2;
        end

    % ===================================================================
    case '6DOF'
    % -------------------------------------------------------------------
    % Full nonlinear 6-DOF dynamics in ECI frame.
    %
    % State: x = [r_ECI(3); v_ECI(3); q_bi(4); omega_b(3)]  (13x1)
    %   r_ECI    [3]  chaser position in ECI [m]
    %   v_ECI    [3]  chaser velocity in ECI [m/s]
    %   q_bi     [4]  scalar-first quaternion: body w.r.t. ICRF
    %   omega_b  [3]  angular velocity in body frame [rad/s]
    %
    % Translational: r'' = -mu/r³ * r + f_body_rotated/m + perturb
    % Rotational (Euler): I*alpha = tau - omega × (I*omega)
    % Kinematic: q_dot = 0.5 * Omega(omega) * q
    % -------------------------------------------------------------------

        if numel(x) ~= 13
            error('6DOF mode requires 13-element state vector');
        end

        r = x(1:3);
        v = x(4:6);
        q = x(7:10);   % scalar-first quaternion
        w = x(11:13);  % body angular rate

        q = q / norm(q);  % normalise

        % Gravitational acceleration (two-body, ECI)
        r_norm = norm(r);
        a_grav = -(mu / r_norm^3) * r;
        aux.acc_grav = a_grav;

        % J2 perturbation (ECI)
        if getparam(params, 'J2_enable', false)
            R_E = 6.371e6;
            J2  = 1.08263e-3;
            z   = r(3);
            fac = (3/2) * J2 * mu * R_E^2 / r_norm^5;
            a_J2 = fac * [(5*z^2/r_norm^2 - 1)*r(1); ...
                          (5*z^2/r_norm^2 - 1)*r(2); ...
                          (5*z^2/r_norm^2 - 3)*r(3)];
            aux.acc_J2 = a_J2;
        else
            a_J2 = zeros(3,1);
        end

        % Atmospheric drag perturbation (ECI, simplified)
        if getparam(params, 'drag_enable', false)
            CD  = getparam(params, 'CD',      2.2);
            A   = getparam(params, 'A_m2',    10);   % m²
            rho = getparam(params, 'rho_kgm3',1e-12); % kg/m³ at ISS altitude
            v_rel = v;  % Ignore Earth rotation (small effect)
            a_drag = -0.5 * CD * A * rho / mass * norm(v_rel) * v_rel;
            aux.acc_drag = a_drag;
        else
            a_drag = zeros(3,1);
        end

        % Rotate applied body force to ECI
        DCM_b2i = quat2dcm_sf(q)';  % body -> ECI  (transpose of ECI->body)
        a_thrust = DCM_b2i * u_force / mass;

        % Translational equations of motion
        xdot(1:3) = v;
        xdot(4:6) = a_grav + a_J2 + a_drag + a_thrust;

        % Attitude kinematics: q_dot = 0.5 * Xi(q) * omega_body
        % q_dot = 0.5 * [  -q_vec'  ] * omega
        %                 [q0*I + q_vec×]
        q0 = q(1); qv = q(2:4);
        Xi = [-qv'; q0*eye(3) + skew(qv)];  % 4x3 matrix
        xdot(7:10) = 0.5 * Xi * w;

        % Rotational dynamics: Euler's equation
        %   I*omega_dot = tau - omega × (I*omega)
        I_b = getparam(params, 'Inertia_kgm2', ...
                        diag([5000 4000 3000]));  % kg⋅m² (typical chaser)

        Iw  = I_b * w;
        xdot(11:13) = I_b \ (u_torque - cross(w, Iw));

    % ===================================================================
    case 'keplerian'
    % -------------------------------------------------------------------
    % Restricted two-body (Keplerian) absolute dynamics
    % State: x = [r_ECI(3); v_ECI(3)]  (6x1)
    % Used for: far-range phasing, initialisation of CW/6DOF
    % -------------------------------------------------------------------

        if numel(x) ~= 6
            error('keplerian mode requires 6-element state vector');
        end

        r = x(1:3);
        v = x(4:6);

        r_norm = norm(r);
        a_grav = -(mu / r_norm^3) * r;
        aux.acc_grav = a_grav;

        % J2 perturbation
        if getparam(params, 'J2_enable', false)
            R_E = 6.371e6;
            J2  = 1.08263e-3;
            z   = r(3);
            fac = (3/2) * J2 * mu * R_E^2 / r_norm^5;
            a_J2 = fac * [(5*z^2/r_norm^2 - 1)*r(1); ...
                          (5*z^2/r_norm^2 - 1)*r(2); ...
                          (5*z^2/r_norm^2 - 3)*r(3)];
        else
            a_J2 = zeros(3,1);
        end

        a_thrust = u_force / mass;

        xdot(1:3) = v;
        xdot(4:6) = a_grav + a_J2 + a_thrust;

    otherwise
        error('RD_RelativeDynamics: Unknown mode ''%s''', mode);
end

end

%% -----------------------------------------------------------------------
%  Coordinate conversion utilities
% -----------------------------------------------------------------------

function S = skew(v)
%SKEW  Returns 3x3 skew-symmetric matrix for cross product: skew(a)*b = a×b
    S = [ 0,    -v(3),  v(2); ...
          v(3),  0,    -v(1); ...
         -v(2),  v(1),  0   ];
end

function DCM = quat2dcm_sf(q)
%QUAT2DCM_SF  Rotation matrix (ECI->body) from scalar-first quaternion
    w = q(1); x = q(2); y = q(3); z = q(4);
    DCM = [1-2*(y^2+z^2),  2*(x*y+w*z),   2*(x*z-w*y); ...
           2*(x*y-w*z),  1-2*(x^2+z^2),   2*(y*z+w*x); ...
           2*(x*z+w*y),   2*(y*z-w*x),  1-2*(x^2+y^2)];
end

function val = getparam(params, field, default)
    if isfield(params, field)
        val = params.(field);
    else
        val = default;
    end
end

%% -----------------------------------------------------------------------
%  CW Analytical Solution (for reference / maneuver planning)
%  Fehse Eq. 3.42 – state transition matrix
% -----------------------------------------------------------------------
function Phi = RD_CW_STM(n, dt)
%RD_CW_STM  Clohessy-Wiltshire State Transition Matrix over time dt
%   Phi = RD_CW_STM(n, dt) returns the 6x6 STM for mean motion n [rad/s]
%   propagating a CW state x=[dx;dy;dz;dxdot;dydot;dzdot] over dt [s].
%
%   Fehse (2003), Eq. 3.42

    s = sin(n*dt);
    c = cos(n*dt);

    Phi = [4-3*c,    0, 0, s/n,        2*(1-c)/n,  0;    ...
           6*(s-n*dt), 1, 0, -2*(1-c)/n, (4*s-3*n*dt)/n, 0; ...
           0,        0, c, 0,           0,          s/n;  ...
           3*n*s,    0, 0, c,           2*s,        0;    ...
          -6*n*(1-c),0, 0, -2*s,        4*c-3,      0;    ...
           0,        0, -n*s, 0,        0,          c];
end
