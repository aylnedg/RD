function [F_cmd, T_cmd, gd] = RD_GuidanceControl(phase_id, nav, target, params)
%RD_GUIDANCECONTROL  Phase-specific guidance and control laws for R&D
%   [F_cmd, T_cmd, gd] = RD_GuidanceControl(phase_id, nav, target, params)
%
%   Computes commanded force (N, body frame) and torque (Nm, body frame)
%   for each mission phase according to Fehse's strategies.
%
%   Guidance strategies (Fehse Ch.5,6,7):
%   Phase 0 – Launch       : Ascent guidance (simplified PD)
%   Phase 1 – Phasing      : Hohmann transfer (impulsive dV commands)
%   Phase 2 – Far R&V      : Radial approach using CW / Lambert
%   Phase 3 – Homing       : V-bar or R-bar approach using CW targeting
%   Phase 4 – Closing      : Forced-motion along approach axis
%   Phase 5 – Final Appr.  : 6-DOF corridor-constrained approach
%   Phase 6 – Docking      : Null-rates + capture engagement
%
%   Inputs:
%     phase_id  integer   Current phase (0–6)
%     nav       struct    Navigation outputs from RD_NavigationFilter:
%                         .rel_pos_LVLH [3x1] m
%                         .rel_vel_LVLH [3x1] m/s
%                         .att_quat     [4x1]
%                         .ang_rate     [3x1] rad/s
%                         .range_m      scalar
%     target    struct    Target state (for outer-loop planning):
%                         .pos_ECI [3x1], .vel_ECI [3x1]
%     params    struct    Guidance/control parameters (overrides defaults)
%
%   Outputs:
%     F_cmd  [3x1] N     Commanded force (LVLH frame for CW; body for 6DOF)
%     T_cmd  [3x1] Nm    Commanded torque (body frame)
%     gd     struct      Guidance diagnostics (ref trajectory, errors, etc.)
%
% Reference:
%   Fehse, W. (2003). Automated Rendezvous and Docking of Spacecraft.
%   Cambridge University Press. Ch. 5,6,7.
%
% See also: RD_NavigationFilter, RD_ActuatorModels, RD_RunAllPhases

% Copyright 2024 - Extended R&D Simulation

F_cmd = zeros(3,1);
T_cmd = zeros(3,1);
gd    = struct('phase', phase_id, 'mode', '', 'error', zeros(3,1), ...
               'ref_pos', zeros(3,1), 'ref_vel', zeros(3,1));

n    = getparam(params, 'n',       1.1e-3);    % mean motion [rad/s]
mass = getparam(params, 'mass_kg', 8000);      % chaser mass [kg]

rel_pos = nav.rel_pos_LVLH;
rel_vel = nav.rel_vel_LVLH;

switch phase_id

    % ===================================================================
    case 0  % Launch & Ascent
    % -------------------------------------------------------------------
    % Simplified ascent guidance: pitch program + thrust vectoring
    % Real implementation uses explicit guidance (e.g., Powered-Explicit)
    % -------------------------------------------------------------------
        gd.mode = 'ascent_pitch_program';

        % Target: circular orbit at desired altitude
        h_target = getparam(params, 'target_alt_m', 400e3); % 400 km
        R_E      = 6.371e6;
        r_target = R_E + h_target;
        mu_E     = 3.986004418e14;
        v_circ   = sqrt(mu_E / r_target);

        % Simple PD: error in altitude and velocity
        if isfield(nav, 'pos_ECI') && isfield(nav, 'vel_ECI')
            r_now = norm(nav.pos_ECI);
            alt_err = r_target - r_now;
            F_cmd = [mass * 0.1 * alt_err; 0; 0];  % Simplified radial thrust
        end

        % Attitude: pitch for trajectory
        T_cmd = zeros(3,1);

    % ===================================================================
    case 1  % Phasing
    % -------------------------------------------------------------------
    % Two-impulse Hohmann-type phasing (Fehse Ch.5.3)
    % Calculates delta-V for each burn based on required phase angle change
    % -------------------------------------------------------------------
        gd.mode = 'phasing_Hohmann';

        dphase = getparam(params, 'phase_angle_deg', 0);  % degrees ahead/behind
        T_orbit = 2*pi/n;

        % Required phasing orbit period for N revolutions
        N_rev  = getparam(params, 'phasing_revs', 2);
        T_phase = T_orbit - dphase*(pi/180)/n/N_rev;

        % Semi-major axis of phasing orbit (from Kepler's 3rd law)
        mu_E   = 3.986004418e14;
        a_phase = (mu_E * (T_phase/(2*pi))^2)^(1/3);
        a_tgt  = mu_E^(1/3) * (T_orbit/(2*pi))^(2/3);

        % Delta-V for Hohmann: at apogee/perigee transitions
        v_circ  = sqrt(mu_E / a_tgt);
        v_phase_p = sqrt(mu_E * (2/a_tgt - 1/a_phase));  % at periapsis
        dV1 = v_phase_p - v_circ;     % First burn (along-track)
        dV2 = v_circ - v_phase_p;     % Second burn (re-circularise)

        % Output as along-track (y-axis LVLH) force command
        % In real mission: output as event-triggered delta-V pulses
        F_cmd = [0; dV1 * mass / 100; 0];  % spread over 100 s burn

        gd.dV1_ms = dV1;
        gd.dV2_ms = dV2;

    % ===================================================================
    case 2  % Far-Range Rendezvous
    % -------------------------------------------------------------------
    % CW-based Lambert guidance for approach to S2 (50 km → 5 km)
    % Uses 2-impulse transfer targeting final position in LVLH
    % Fehse Ch.5.5 – radial approach corridor
    % -------------------------------------------------------------------
        gd.mode = 'far_range_Lambert_CW';

        % Target: arrive at waypoint S2 (~5 km behind, on V-bar)
        r_target_LVLH = getparam(params, 'wp_S2_m', [-5000; 0; 0]);  % radial
        t_transfer    = getparam(params, 'transfer_time_s', 3600);

        % CW two-impulse transfer: solve for initial velocity
        [dV1, dV2] = CW_lambert_2impulse(rel_pos, r_target_LVLH, n, t_transfer);

        F_cmd = mass * dV1 / max(getparam(params,'burn_duration_s',60), 1);
        gd.dV1_ms = dV1;
        gd.dV2_ms = dV2;

    % ===================================================================
    case 3  % Homing / Close-Range Rendezvous
    % -------------------------------------------------------------------
    % V-bar or R-bar approach using CW-based guidance
    % Fehse Ch.6.3 – forced motion approach trajectories
    % -------------------------------------------------------------------
        approach_axis = getparam(params, 'approach_axis', 'Vbar');
        gd.mode = ['homing_' approach_axis];

        % Final waypoint S3 (~200 m from target)
        if strcmp(approach_axis, 'Vbar')
            r_wp = [0; -200; 0];   % 200 m behind on V-bar
        else  % R-bar
            r_wp = [-200; 0; 0];   % 200 m below on R-bar
        end
        r_wp = getparam(params, 'wp_S3_m', r_wp);

        % PD guidance law with CW feed-forward
        Kp = getparam(params, 'Kp_homing', 0.002);  % pos gain
        Kd = getparam(params, 'Kd_homing', 0.05);   % vel gain

        pos_err = r_wp - rel_pos;
        vel_ref  = Kp * pos_err;           % Proportional velocity reference
        vel_err  = vel_ref - rel_vel;

        % CW gravity gradient compensation (Fehse Eq. 3.38 inverse)
        cw_ff = [3*n^2*rel_pos(1) + 2*n*rel_vel(2); ...
                 -2*n*rel_vel(1);                    ...
                 -n^2*rel_pos(3)];

        a_cmd = Kd * vel_err - cw_ff;  % Note: subtract to cancel CW terms
        F_cmd = mass * a_cmd;

        gd.error   = pos_err;
        gd.ref_pos = r_wp;
        gd.ref_vel = vel_ref;

    % ===================================================================
    case 4  % Closing / Forced-Motion Approach
    % -------------------------------------------------------------------
    % Forced motion along V-bar (or R-bar) with constant velocity
    % Safety: corridor constraints on lateral deviation (Fehse Ch.6.4)
    % -------------------------------------------------------------------
        approach_axis = getparam(params, 'approach_axis', 'Vbar');
        gd.mode = ['closing_forced_' approach_axis];

        v_approach = getparam(params, 'approach_speed_ms', -0.05);  % m/s (-y = closing)
        corridor_r = getparam(params, 'corridor_radius_m', 2);       % m

        % Reference: constant-velocity closing on approach axis
        if strcmp(approach_axis, 'Vbar')
            v_ref = [0; v_approach; 0];
            pos_ref = [0; rel_pos(2); 0];  % Stay on V-bar
        else  % R-bar
            v_ref   = [v_approach; 0; 0];
            pos_ref = [rel_pos(1); 0; 0];
        end

        % Lateral deviation controller (corridor constraint)
        if strcmp(approach_axis, 'Vbar')
            lat_dev   = [rel_pos(1); rel_pos(3)];      % radial + cross-track
            lat_vel   = [rel_vel(1); rel_vel(3)];
        else
            lat_dev   = [rel_pos(2); rel_pos(3)];
            lat_vel   = [rel_vel(2); rel_vel(3)];
        end

        Kp_lat = getparam(params, 'Kp_lat', 0.01);
        Kd_lat = getparam(params, 'Kd_lat', 0.2);
        a_lat  = -Kp_lat * lat_dev - Kd_lat * lat_vel;

        % Along-axis velocity control
        Kd_ax = getparam(params, 'Kd_ax', 0.5);
        if strcmp(approach_axis, 'Vbar')
            a_ax = Kd_ax * (v_approach - rel_vel(2));
            a_cmd = [a_lat(1); a_ax; a_lat(2)];
        else
            a_ax = Kd_ax * (v_approach - rel_vel(1));
            a_cmd = [a_ax; a_lat(1); a_lat(2)];
        end

        % CW feed-forward
        cw_ff = [3*n^2*rel_pos(1) + 2*n*rel_vel(2); ...
                 -2*n*rel_vel(1); ...
                 -n^2*rel_pos(3)];
        F_cmd = mass * (a_cmd - cw_ff);

        % Corridor violation warning
        if norm(lat_dev) > corridor_r
            gd.corridor_violation = true;
        else
            gd.corridor_violation = false;
        end

        gd.error   = rel_pos - pos_ref;
        gd.ref_vel = v_ref;

    % ===================================================================
    case 5  % Final Approach
    % -------------------------------------------------------------------
    % 6-DOF approach: position + attitude control
    % Corridor constraint: ±15 deg cone around docking axis
    % Fehse Ch.6.5 – docking port alignment
    % -------------------------------------------------------------------
        gd.mode = 'final_approach_6DOF';

        v_approach = getparam(params, 'approach_speed_ms', -0.01); % m/s
        Kp_pos = getparam(params, 'Kp_pos', 0.02);
        Kd_pos = getparam(params, 'Kd_pos', 0.3);
        Kp_att = getparam(params, 'Kp_att', 0.1);
        Kd_att = getparam(params, 'Kd_att', 0.5);

        % Position control (approach axis = -y LVLH for V-bar)
        v_ref = [0; v_approach; 0];
        lat_err = [rel_pos(1); rel_pos(3)];
        lat_vel = [rel_vel(1); rel_vel(3)];

        a_lat  = -Kp_pos * lat_err - Kd_pos * lat_vel;
        a_ax   = Kd_pos * (v_approach - rel_vel(2));
        a_cmd  = [a_lat(1); a_ax; a_lat(2)];

        % CW feed-forward
        cw_ff = [3*n^2*rel_pos(1) + 2*n*rel_vel(2); ...
                 -2*n*rel_vel(1); ...
                 -n^2*rel_pos(3)];
        F_cmd = mass * (a_cmd - cw_ff);

        % Attitude control: align docking port with target port
        % Target: chaser +y body axis aligned with docking direction (-y LVLH)
        q_ref = getparam(params, 'q_dock_ref', [1;0;0;0]);  % Desired attitude

        if isfield(nav, 'att_quat')
            q_err = quatmultiply_sf(q_ref, quatconj_sf(nav.att_quat));
            att_err_vec = q_err(2:4) * sign(q_err(1));  % Error rotation vector
        else
            att_err_vec = zeros(3,1);
        end

        T_cmd = Kp_att * att_err_vec - Kd_att * nav.ang_rate;

        % Cone corridor check (15 deg half-angle)
        if norm(lat_err) / max(abs(rel_pos(2)), 1e-3) > tand(15)
            gd.corridor_violation = true;
        else
            gd.corridor_violation = false;
        end

        gd.error   = rel_pos;
        gd.ref_vel = v_ref;

    % ===================================================================
    case 6  % Docking
    % -------------------------------------------------------------------
    % Null-velocity control + capture mechanism engagement
    % Fehse Ch.7 – contact dynamics and rigidisation
    % -------------------------------------------------------------------
        gd.mode = 'docking_null_rates';

        Kd_null = getparam(params, 'Kd_null', 0.5);
        Kd_att  = getparam(params, 'Kd_att',  0.5);

        % Null any residual translational velocity
        F_cmd = -Kd_null * mass * rel_vel;

        % Null any angular rates
        if isfield(nav, 'ang_rate')
            T_cmd = -Kd_att * nav.ang_rate;
        end

        gd.ref_pos = zeros(3,1);
        gd.ref_vel = zeros(3,1);

    otherwise
        warning('RD_GuidanceControl: Unknown phase_id %d', phase_id);
end

end

%% -----------------------------------------------------------------------
%  CW Lambert Two-Impulse Solver
%  Fehse Appendix A1 – CW transfer between two points in time dt
% -----------------------------------------------------------------------
function [dV1, dV2] = CW_lambert_2impulse(r0, rf, n, dt)
%CW_LAMBERT_2IMPULSE  Two-impulse transfer in CW frame
%   Solves for initial velocity required to transfer from r0 to rf in time dt

    s = sin(n*dt); c = cos(n*dt);

    % CW state transition matrix (position partitions)
    % x(t) = Phi11*x0 + Phi12*v0  →  v0 = Phi12^{-1}*(rf - Phi11*r0)
    Phi11 = [4-3*c,      0, 0; ...
             6*(s-n*dt), 1, 0; ...
             0,          0, c];

    Phi12 = [s/n,         2*(1-c)/n,       0; ...
             -2*(1-c)/n, (4*s-3*n*dt)/n,  0; ...
             0,           0,              s/n];

    Phi21 = [3*n*s,      0, 0; ...
             -6*n*(1-c), 0, 0; ...
             0,          0,-n*s];

    Phi22 = [c,  2*s,     0; ...
             -2*s, 4*c-3, 0; ...
             0,   0,      c];

    % Solve for required initial velocity
    rhs = rf - Phi11 * r0;
    cond_check = cond(Phi12);
    if cond_check > 1e10
        warning('CW_lambert: ill-conditioned matrix (cond=%.2e)', cond_check);
        v0_required = pinv(Phi12) * rhs;
    else
        v0_required = Phi12 \ rhs;
    end

    dV1 = v0_required;  % Required velocity change at r0

    % Final velocity (at rf)
    vf = Phi21 * r0 + Phi22 * v0_required;
    dV2 = -vf;  % Kill residual velocity at rf (for rendezvous)
end

%% -----------------------------------------------------------------------
%  Quaternion helpers
% -----------------------------------------------------------------------

function q_ab = quatmultiply_sf(q_a, q_b)
%QUATMULTIPLY_SF  Quaternion product, scalar-first
    w1 = q_a(1); v1 = q_a(2:4);
    w2 = q_b(1); v2 = q_b(2:4);
    q_ab = [w1*w2 - dot(v1,v2); w1*v2 + w2*v1 + cross(v1,v2)];
end

function q_c = quatconj_sf(q)
    q_c = [q(1); -q(2:4)];
end

function val = getparam(params, field, default)
    if isfield(params, field)
        val = params.(field);
    else
        val = default;
    end
end
