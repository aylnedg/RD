function out = RD_ActuatorModels(actuatorType, cmd, state, params)
%RD_ACTUATORMODELS  Actuator force/torque models for all R&D phases
%   out = RD_ActuatorModels(actuatorType, cmd, state, params)
%
%   ACTUATOR TYPES (actuatorType string):
%     'MainEngine'           - High-thrust bipropellant main engine
%     'RCS_Thrusters'        - Reaction Control System thruster set (8–16 thrusters)
%     'RCS_Thrusters_Fine'   - Fine RCS for final approach (smaller thrust)
%     'ReactionWheels'       - Momentum-exchange attitude actuators
%     'MagneticTorquer'      - Magnetic torque rods (magnetorquers)
%     'DockingMechanism'     - Capture ring / drogue mechanism
%
%   cmd    : Command struct. Fields depend on actuator (see below).
%   state  : Vehicle state struct (att_quat, ang_rate, pos_ECI, etc.)
%   params : (optional) struct to override default parameters.
%
%   OUTPUTS (out struct):
%     .force_body_N   [3x1] N   Net force in body frame
%     .torque_body_Nm [3x1] Nm  Net torque about CoM in body frame
%     .mass_flow_kgs  scalar    Propellant mass flow rate (kg/s)
%     .power_W        scalar    Electrical power consumed (W)
%     .status         string    'nominal' | 'saturated' | 'off'
%
% Reference:
%   Fehse, W. (2003). Automated Rendezvous and Docking of Spacecraft.
%   Cambridge University Press. ISBN 0-521-82492-3.
%
% See also: RD_Fehse_PhaseConfig, RD_SensorModels, RD_GuidanceControl

% Copyright 2024 - Extended R&D Simulation

if nargin < 4 || isempty(params)
    params = struct();
end

% Initialise output
out.force_body_N   = zeros(3,1);
out.torque_body_Nm = zeros(3,1);
out.mass_flow_kgs  = 0;
out.power_W        = 0;
out.status         = 'nominal';

switch actuatorType

    % ===================================================================
    case 'MainEngine'
    % -------------------------------------------------------------------
    % High-thrust bipropellant main engine
    % Typical: 490 N bi-prop (e.g., ArianeGroup S10-18)
    % Isp: 321 s (NTO/MMH)
    % Used for: large delta-V maneuvers in phasing, far-range rendezvous
    % Command: cmd.dV_ECI [3x1] m/s (desired delta-V in ECI)
    %          cmd.burn_duration_s  scalar
    % -------------------------------------------------------------------
        F_max  = getparam(params, 'max_thrust_N', 490);   % N
        Isp    = getparam(params, 'Isp_s',        321);   % s
        g0     = 9.80665;  % m/s²

        if isfield(cmd, 'dV_ECI') && norm(cmd.dV_ECI) > 1e-6
            % Convert desired delta-V direction to body frame
            DCM_b2e = quat2dcm_sf(state.att_quat)';  % body to ECI
            dv_body = DCM_b2e \ cmd.dV_ECI;           % ECI -> body (approx)
            thrust_dir = dv_body / norm(dv_body);

            % Commanded thrust (assume full thrust burn)
            thrust_body = F_max * thrust_dir;

            % Apply misalignment (0.1 deg typical)
            sigma_mis = getparam(params, 'misalign_rad', 1.745e-3);
            noise_perp = sigma_mis * randn(3,1);
            noise_perp(3) = -dot(thrust_dir(1:2), noise_perp(1:2)) / thrust_dir(3);
            thrust_body = thrust_body + F_max * noise_perp;

            out.force_body_N   = thrust_body;
            out.torque_body_Nm = zeros(3,1);  % Assume CoT at CoM
            out.mass_flow_kgs  = F_max / (Isp * g0);
            out.power_W        = 0;  % Pressure-fed, no extra power
            out.status         = 'nominal';
        else
            out.status = 'off';
        end

    % ===================================================================
    case 'RCS_Thrusters'
    % -------------------------------------------------------------------
    % Reaction Control System – 6-DOF thruster set
    % Configuration: 12 thrusters in 2 clusters (+x/-x each with ±y,±z pairs)
    % Typical: 22 N bi-prop (NTO/MMH), Isp ~220 s
    % Used for: attitude control + small orbital maneuvers (all close phases)
    % Command: cmd.force_body_N  [3x1] desired force in body frame
    %          cmd.torque_body_Nm [3x1] desired torque
    % -------------------------------------------------------------------
        F_max     = getparam(params, 'max_thrust_N', 22);   % N per thruster
        Isp       = getparam(params, 'Isp_s',        220);  % s
        n_thrusters = getparam(params, 'n_thrusters', 12);
        g0        = 9.80665;

        % Pseudo-inverse allocation (simplified – axis-aligned thrusters)
        % Each axis can provide F_max in both directions
        F_max_total = F_max * floor(n_thrusters / 6);  % thrusters per axis-dir

        if isfield(cmd, 'force_body_N')
            F_cmd = cmd.force_body_N;
        else
            F_cmd = zeros(3,1);
        end
        if isfield(cmd, 'torque_body_Nm')
            T_cmd = cmd.torque_body_Nm;
        else
            T_cmd = zeros(3,1);
        end

        % Saturate force and torque
        sat_F = min(1, F_max_total ./ max(abs(F_cmd), 1e-9));
        sat_T_max = F_max * 0.5;  % 0.5 m moment arm
        sat_T = min(1, sat_T_max ./ max(abs(T_cmd), 1e-9));
        sat = min(min(sat_F), min(sat_T));

        out.force_body_N   = F_cmd * sat;
        out.torque_body_Nm = T_cmd * sat;

        % Mass flow (proportional to total impulse)
        total_thrust = norm(F_cmd) + norm(T_cmd)/0.5;
        out.mass_flow_kgs  = total_thrust / (Isp * g0);
        out.power_W        = 50;  % valve driver power

        if sat < 0.99
            out.status = 'saturated';
        end

    % ===================================================================
    case 'RCS_Thrusters_Fine'
    % -------------------------------------------------------------------
    % Fine RCS for final approach: smaller thrusters for precision
    % Typical: 4 N monoprop (Hydrazine), Isp ~220 s
    % Used for: precise 6-DOF control in phases 5,6
    % Command: same as RCS_Thrusters
    % -------------------------------------------------------------------
        F_max     = getparam(params, 'max_thrust_N', 4);    % N per thruster (small)
        Isp       = getparam(params, 'Isp_s',        220);
        n_thrusters = getparam(params, 'n_thrusters', 12);
        g0        = 9.80665;

        F_max_total = F_max * floor(n_thrusters / 6);

        F_cmd = isfield_or_zero(cmd, 'force_body_N',   3);
        T_cmd = isfield_or_zero(cmd, 'torque_body_Nm', 3);

        sat_F = min(1, F_max_total ./ max(abs(F_cmd), 1e-9));
        sat_T = min(1, (F_max * 0.3) ./ max(abs(T_cmd), 1e-9));
        sat   = min(min(sat_F), min(sat_T));

        out.force_body_N   = F_cmd * sat;
        out.torque_body_Nm = T_cmd * sat;
        out.mass_flow_kgs  = (norm(F_cmd) + norm(T_cmd)/0.3) / (Isp * g0);
        out.power_W        = 20;

        if sat < 0.99
            out.status = 'saturated';
        end

    % ===================================================================
    case 'ReactionWheels'
    % -------------------------------------------------------------------
    % Reaction wheel assembly (RWA) – 4-wheel pyramid configuration
    % Typical: Honeywell HR14-100 / Airbus W18
    % Max torque: 0.2 Nm per wheel, Max momentum: 100 Nms
    % Used for: precise attitude control in all phases
    % Command: cmd.torque_body_Nm [3x1] desired torque in body frame
    %          cmd.wheel_speed_rads [4x1] current wheel speeds
    % -------------------------------------------------------------------
        T_max_wheel = getparam(params, 'max_torque_Nm_wheel', 0.2);   % Nm
        H_max_wheel = getparam(params, 'max_momentum_Nms',    100);   % Nms
        power_coeff = getparam(params, 'power_W_per_Nm',      10);    % W/Nm

        % 4-wheel pyramid mounting matrix [A] (torque distribution)
        % Wheels tilted 54.74 deg from spin axis (tetrahedral config)
        beta = 54.74 * pi/180;
        A = [cos(beta),  cos(beta),  cos(beta),  cos(beta); ...
             -sin(beta),  sin(beta), -sin(beta),  sin(beta); ...
              0,          0,          0,          0];
        A(3,:) = cross(A(1,:), A(2,:));  % Ensure orthogonality (approx)
        A = [sin(0)*cos(0), sin(0)*cos(pi/2), sin(0)*cos(pi), sin(0)*cos(3*pi/2); ...
             sin(0)*sin(0), sin(0)*sin(pi/2), sin(0)*sin(pi), sin(0)*sin(3*pi/2); ...
             cos(beta),     cos(beta),         cos(beta),      cos(beta)];
        % Simplified: direct 3-axis control with saturation
        T_cmd = isfield_or_zero(cmd, 'torque_body_Nm', 3);

        T_max_body = T_max_wheel * 3;  % 3 wheels contribute per axis (approx)
        sat = min(1, T_max_body ./ max(abs(T_cmd), 1e-9));
        sat_scalar = min(sat);

        out.force_body_N   = zeros(3,1);
        out.torque_body_Nm = T_cmd * sat_scalar;
        out.mass_flow_kgs  = 0;  % Momentum exchange, no propellant
        out.power_W        = power_coeff * norm(T_cmd) + 5;  % 5 W idle

        if sat_scalar < 0.99
            out.status = 'saturated';
        end

        % Momentum saturation check
        if isfield(cmd, 'wheel_speed_rads')
            I_wheel = 0.05;  % kg⋅m² typical wheel inertia
            H_wheels = I_wheel * cmd.wheel_speed_rads;
            if any(abs(H_wheels) > H_max_wheel * 0.9)
                out.status = 'momentum_saturated';
                % Desaturation needed via magnetorquer or RCS
            end
        end

    % ===================================================================
    case 'MagneticTorquer'
    % -------------------------------------------------------------------
    % Magnetic torque rods (magnetorquers) – 3-axis set
    % Typical: ZARM MTQ-800 (800 Am² dipole)
    % Torque: T = m × B (depends on magnetic field)
    % Used for: attitude control + wheel desaturation in LEO phases 0,1,2
    % Command: cmd.dipole_Am2 [3x1] magnetic dipole moment vector
    %          state.mag_vec_body [3x1] T  magnetic field in body frame
    % -------------------------------------------------------------------
        m_max = getparam(params, 'max_dipole_Am2', 800);  % Am²

        m_cmd = isfield_or_zero(cmd, 'dipole_Am2', 3);

        % Saturate dipole
        if norm(m_cmd) > m_max
            m_cmd = m_cmd * (m_max / norm(m_cmd));
            out.status = 'saturated';
        end

        % Torque = m × B (cross product with local field)
        if isfield(state, 'mag_vec_body')
            B = state.mag_vec_body;
        else
            B = [2e-5; 1e-5; 4e-5];  % Typical LEO B-field [T]
        end

        out.force_body_N   = zeros(3,1);
        out.torque_body_Nm = cross(m_cmd, B);
        out.mass_flow_kgs  = 0;
        out.power_W        = 10 * norm(m_cmd) / m_max;  % 10 W max power

    % ===================================================================
    case 'DockingMechanism'
    % -------------------------------------------------------------------
    % Active docking mechanism: capture ring / probe-and-drogue
    % Models capture latch, attenuation springs, and rigidisation locks
    % Used in: Phase 6 (docking)
    % Command: cmd.capture_enable  logical
    %          cmd.rigidise_enable logical
    %          state.rel_pos_LVLH [3x1] m
    %          state.rel_vel_LVLH [3x1] m/s
    % -------------------------------------------------------------------
        k_spring  = getparam(params, 'spring_stiffness_Nm',  200);   % N/m
        c_damper  = getparam(params, 'damper_coeff_Nms',      50);   % Ns/m
        capture_r = getparam(params, 'capture_radius_m',       0.15); % m

        rel_pos = isfield_or_zero(cmd, 'rel_pos_m',   3);
        rel_vel = isfield_or_zero(cmd, 'rel_vel_ms',  3);

        if isfield(state, 'rel_pos_LVLH')
            rel_pos = state.rel_pos_LVLH;
        end
        if isfield(state, 'rel_vel_LVLH')
            rel_vel = state.rel_vel_LVLH;
        end

        dist = norm(rel_pos);

        if isfield(cmd, 'capture_enable') && cmd.capture_enable && ...
                dist < capture_r
            % Spring-damper contact model
            compression = max(0, capture_r - dist);
            F_spring = k_spring * compression;
            F_damper = c_damper * dot(-rel_vel, rel_pos/max(dist,1e-3));
            F_total  = (F_spring + F_damper) * rel_pos / max(dist, 1e-3);

            out.force_body_N   = -F_total;  % Reaction on chaser
            out.torque_body_Nm = zeros(3,1);
            out.status         = 'capture_active';
        elseif isfield(cmd, 'rigidise_enable') && cmd.rigidise_enable
            out.status = 'rigidised';
        else
            out.status = 'standby';
        end

        out.mass_flow_kgs = 0;
        out.power_W       = 100;  % Mechanism drive power

    otherwise
        error('RD_ActuatorModels: Unknown actuator type ''%s''', actuatorType);
end

end

%% -----------------------------------------------------------------------
%  Local helpers
% -----------------------------------------------------------------------

function val = getparam(params, field, default)
    if isfield(params, field)
        val = params.(field);
    else
        val = default;
    end
end

function v = isfield_or_zero(s, field, n)
    if isfield(s, field)
        v = s.(field);
    else
        v = zeros(n, 1);
    end
end

function DCM = quat2dcm_sf(q)
%QUAT2DCM_SF  Rotation matrix from scalar-first quaternion (body <- ECI)
    w = q(1); x = q(2); y = q(3); z = q(4);
    DCM = [1-2*(y^2+z^2),   2*(x*y+w*z),   2*(x*z-w*y); ...
           2*(x*y-w*z),   1-2*(x^2+z^2),   2*(y*z+w*x); ...
           2*(x*z+w*y),   2*(y*z-w*x),   1-2*(x^2+y^2)];
end
