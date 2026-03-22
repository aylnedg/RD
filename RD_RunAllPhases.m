%RD_RUNALLPHASES  Master simulation runner for all Fehse R&D phases
%
%   Runs a pure-MATLAB (no Simulink GUI required) Monte Carlo-ready
%   simulation of the complete chaser spacecraft mission from launch
%   to docking, cycling through all 7 phases defined by Fehse (2003).
%
%   Outputs:
%     sim_log  struct array  One entry per simulated time step containing:
%                        .t, .phase_id, .phase_name,
%                        .rel_pos_LVLH, .rel_vel_LVLH,
%                        .abs_pos_ECI, .abs_vel_ECI,
%                        .att_quat, .ang_rate,
%                        .range_m, .F_cmd, .T_cmd,
%                        .sensor_validity  (struct of booleans per sensor)
%                        .propellant_kg    (accumulated)
%
%   Usage:
%     sim_log = RD_RunAllPhases();
%     RD_PlotResults(sim_log);
%
% Reference:
%   Fehse, W. (2003). Automated Rendezvous and Docking of Spacecraft.
%   Cambridge University Press. ISBN 0-521-82492-3.
%
% See also: RD_BuildExtendedModel, RD_Fehse_PhaseConfig,
%           RD_SensorModels, RD_ActuatorModels, RD_GuidanceControl,
%           RD_NavigationFilter, RD_RelativeDynamics

% Copyright 2024 - Extended R&D Simulation

clear; clc;
fprintf('==========================================================\n');
fprintf('  Spacecraft Rendezvous & Docking – All Phases Simulation\n');
fprintf('  Based on Fehse (2003) + MathWorks Aerospace Blockset\n');
fprintf('==========================================================\n\n');

%% -----------------------------------------------------------------------
%  1. Simulation parameters
% -----------------------------------------------------------------------
mu_E  = 3.986004418e14;  % m³/s²
R_E   = 6.371e6;          % m
h_ISS = 410e3;            % m  (ISS altitude)
a_ISS = R_E + h_ISS;

n     = sqrt(mu_E / a_ISS^3);   % mean motion ~1.1e-3 rad/s
T_orb = 2*pi / n;               % orbital period ~92 min

% Target initial state (circular orbit, ECI)
r_target_0 = [a_ISS; 0; 0];
v_target_0 = [0; sqrt(mu_E/a_ISS); 0];

% Chaser initial state: 300 km behind (along-track, LVLH), same altitude
r_chaser_0 = r_target_0 - [0; 300e3; 0];  % 300 km behind on V-bar
v_chaser_0 = v_target_0;

% Chaser vehicle properties
mass_0    = 8500;  % kg  (incl. propellant)
I_kgm2    = diag([5000; 4000; 3000]);  % inertia tensor

% Simulation time steps per phase
dt_far    = 60;    % s  Phasing, far R&V (1 min step)
dt_close  = 10;    % s  Homing, closing  (10 s step)
dt_final  = 1;     % s  Final approach, docking (1 s step)

%% -----------------------------------------------------------------------
%  2. Load phase configuration
% -----------------------------------------------------------------------
cfg = RD_Fehse_PhaseConfig();
fprintf('Loaded %d mission phases from Fehse (2003) configuration.\n\n', numel(cfg));
print_phase_summary(cfg);

%% -----------------------------------------------------------------------
%  3. Initial state vectors
% -----------------------------------------------------------------------
% 6-state Keplerian: [r_ECI(3); v_ECI(3)]
x_kep_chaser = [r_chaser_0; v_chaser_0];
x_kep_target = [r_target_0; v_target_0];

% 6-state CW relative: [dx; dy; dz; dxdot; dydot; dzdot] in LVLH
% Initially 300 km behind on V-bar, zero relative velocity
x_cw = [0; -300e3; 0; 0; 0; 0];

% 13-state 6DOF: [r_ECI(3); v_ECI(3); q_bi(4); omega_b(3)]
att0 = [1; 0; 0; 0];  % identity quaternion (body = ICRF)
w0   = zeros(3,1);
x_6dof = [r_chaser_0; v_chaser_0; att0; w0];

% Navigation filter state
x_nav = zeros(12,1);
P_nav = eye(12) * 1e6;

% Propellant bookkeeping
prop_remaining = 2000;  % kg  initial propellant mass

%% -----------------------------------------------------------------------
%  4. Phase-by-phase simulation
% -----------------------------------------------------------------------
sim_log   = struct('t',{},'phase_id',{},'phase_name',{},'rel_pos_LVLH',{}, ...
               'rel_vel_LVLH',{},'abs_pos_ECI',{},'abs_vel_ECI',{}, ...
               'att_quat',{},'ang_rate',{},'range_m',{},'F_cmd',{}, ...
               'T_cmd',{},'sensor_meas',{},'propellant_kg',{});  % time-series sim_log
t_now = 0;    % s

% Sun and magnetic field (simplified: constant direction for this demo)
sun_vec_ECI = normalise([1; 0.1; 0.05]);    % approx sun direction in ECI
mag_vec_ECI = [2e-5; 1e-5; 4e-5];          % T  Earth B-field at ISS orbit

%% ----  PHASE 0: Launch  (t = 0 → 600 s) ------
fprintf('\n--- PHASE 0: Launch & Ascent ---\n');
t_phase = 0 : dt_final : 600;
ph_cfg  = cfg(1);  % phase 0

for i = 1:numel(t_phase)
    t = t_phase(i);
    trueState = build_true_state(x_kep_chaser, x_kep_target, ...
                                  [1;0;0;0], zeros(3,1), zeros(3,1), ...
                                  sun_vec_ECI, mag_vec_ECI, ...
                                  0, false);

    meas = collect_measurements(ph_cfg, trueState);
    nav  = struct('rel_pos_LVLH', zeros(3,1), 'rel_vel_LVLH', zeros(3,1), ...
                  'range_m', norm(x_cw(1:3)), 'att_quat', [1;0;0;0], ...
                  'ang_rate', zeros(3,1));

    [F_cmd, T_cmd, gd] = RD_GuidanceControl(0, nav, struct(), ...
                           struct('n', n, 'mass_kg', mass_0));

    params_dyn = struct('mode','keplerian','mu',mu_E,'mass_kg',mass_0, ...
                        'J2_enable',false,'drag_enable',false);
    xdot = RD_RelativeDynamics(t, x_kep_chaser, F_cmd, T_cmd, params_dyn);
    x_kep_chaser = x_kep_chaser + xdot * dt_final;

    sim_log = append_log(sim_log, t, 0, 'Launch', zeros(3,1), zeros(3,1), ...
                     x_kep_chaser(1:3), x_kep_chaser(4:6), ...
                     [1;0;0;0], zeros(3,1), norm(x_kep_chaser(1:3)-x_kep_target(1:3)), ...
                     F_cmd, T_cmd, meas, prop_remaining);
end
t_now = t_now + 600;
fprintf('  Phase 0 complete. t = %.0f s\n', t_now);

%% ----  PHASE 1: Phasing  (300 km → 50 km) ------
fprintf('\n--- PHASE 1: Phasing (300 km → 50 km) ---\n');
ph_cfg = cfg(2);
t_phase = 0 : dt_far : T_orb * 4;  % up to 4 orbital periods

for i = 1:numel(t_phase)
    t = t_now + t_phase(i);

    % Relative state in LVLH
    rel_pos_LVLH = ECI_to_LVLH(x_kep_chaser(1:3) - x_kep_target(1:3), ...
                                x_kep_target(1:3), x_kep_target(4:6));
    rel_vel_LVLH = ECI_to_LVLH(x_kep_chaser(4:6) - x_kep_target(4:6), ...
                                x_kep_target(1:3), x_kep_target(4:6));
    range = norm(rel_pos_LVLH);

    if range < 50e3, break; end  % Transition to Phase 2

    trueState = build_true_state(x_kep_chaser, x_kep_target, ...
                                  [1;0;0;0], zeros(3,1), [0;range;0], ...
                                  sun_vec_ECI, mag_vec_ECI, ...
                                  0, false);
    trueState.rel_pos_LVLH = rel_pos_LVLH;
    trueState.rel_vel_LVLH = rel_vel_LVLH;

    meas = collect_measurements(ph_cfg, trueState);
    nav  = struct('rel_pos_LVLH', rel_pos_LVLH, 'rel_vel_LVLH', rel_vel_LVLH, ...
                  'range_m', range, 'att_quat', [1;0;0;0], 'ang_rate', zeros(3,1));

    [F_cmd, T_cmd, gd] = RD_GuidanceControl(1, nav, struct(), ...
                           struct('n', n, 'mass_kg', mass_0, ...
                                  'phase_angle_deg', 30, 'phasing_revs', 2));

    % Propagate both spacecraft
    params_dyn = struct('mode','keplerian','mu',mu_E,'mass_kg',mass_0, ...
                        'J2_enable',true,'drag_enable',false);
    x_kep_chaser = rk4(x_kep_chaser, F_cmd, T_cmd, t, dt_far, params_dyn);
    x_kep_target = rk4(x_kep_target, zeros(3,1), zeros(3,1), t, dt_far, ...
                        struct('mode','keplerian','mu',mu_E,'mass_kg',20000,...
                               'J2_enable',true,'drag_enable',false));

    prop_remaining = prop_remaining - norm(F_cmd) * dt_far / (300*9.81);

    sim_log = append_log(sim_log, t, 1, 'Phasing', rel_pos_LVLH, rel_vel_LVLH, ...
                     x_kep_chaser(1:3), x_kep_chaser(4:6), ...
                     [1;0;0;0], zeros(3,1), range, F_cmd, T_cmd, meas, prop_remaining);
end
t_now = t_now + t_phase(min(i, numel(t_phase)));
fprintf('  Phase 1 complete. Range = %.1f km. t = %.0f s\n', range/1e3, t_now);

%% ----  PHASE 2: Far-Range Rendezvous  (50 km → 5 km) ------
fprintf('\n--- PHASE 2: Far-Range Rendezvous (50 km → 5 km) ---\n');
ph_cfg = cfg(3);
x_cw   = [rel_pos_LVLH; rel_vel_LVLH];  % Switch to CW state

for i = 1:10000
    t = t_now + (i-1)*dt_far;
    range = norm(x_cw(1:3));
    if range < 5e3, break; end

    trueState = build_true_state(x_kep_chaser, x_kep_target, ...
                                  [1;0;0;0], zeros(3,1), zeros(3,1), ...
                                  sun_vec_ECI, mag_vec_ECI, ...
                                  0, false);
    trueState.rel_pos_LVLH = x_cw(1:3);
    trueState.rel_vel_LVLH = x_cw(4:6);

    meas = collect_measurements(ph_cfg, trueState);
    nav  = struct('rel_pos_LVLH', x_cw(1:3), 'rel_vel_LVLH', x_cw(4:6), ...
                  'range_m', range, 'att_quat', [1;0;0;0], 'ang_rate', zeros(3,1));

    [F_cmd, T_cmd] = RD_GuidanceControl(2, nav, struct(), ...
                      struct('n', n, 'mass_kg', mass_0, ...
                             'wp_S2_m', [-5000; 0; 0], ...
                             'transfer_time_s', 3600));

    params_cw = struct('mode','CW','n',n,'mass_kg',mass_0,'J2_enable',false);
    xdot_cw = RD_RelativeDynamics(t, x_cw, F_cmd, T_cmd, params_cw);
    x_cw = x_cw + xdot_cw * dt_far;

    prop_remaining = prop_remaining - norm(F_cmd) * dt_far / (220*9.81);

    sim_log = append_log(sim_log, t, 2, 'FarRendezvous', x_cw(1:3), x_cw(4:6), ...
                     x_kep_chaser(1:3), x_kep_chaser(4:6), ...
                     [1;0;0;0], zeros(3,1), range, F_cmd, T_cmd, meas, prop_remaining);
end
t_now = t_now + (i-1)*dt_far;
fprintf('  Phase 2 complete. Range = %.0f m. t = %.0f s\n', range, t_now);

%% ----  PHASE 3: Homing  (5 km → 200 m) ------
fprintf('\n--- PHASE 3: Homing (5 km → 200 m) ---\n');
ph_cfg = cfg(4);

for i = 1:50000
    t = t_now + (i-1)*dt_close;
    range = norm(x_cw(1:3));
    if range < 200, break; end

    trueState = build_true_state(x_kep_chaser, x_kep_target, ...
                                  [1;0;0;0], zeros(3,1), zeros(3,1), ...
                                  sun_vec_ECI, mag_vec_ECI, ...
                                  0, false);
    trueState.rel_pos_LVLH = x_cw(1:3);
    trueState.rel_vel_LVLH = x_cw(4:6);

    meas = collect_measurements(ph_cfg, trueState);
    nav  = struct('rel_pos_LVLH', x_cw(1:3), 'rel_vel_LVLH', x_cw(4:6), ...
                  'range_m', range, 'att_quat', [1;0;0;0], 'ang_rate', zeros(3,1));

    [F_cmd, T_cmd] = RD_GuidanceControl(3, nav, struct(), ...
                      struct('n', n, 'mass_kg', mass_0, ...
                             'approach_axis', 'Vbar', ...
                             'wp_S3_m', [0; -200; 0]));

    params_cw = struct('mode','CW','n',n,'mass_kg',mass_0,'J2_enable',false);
    xdot_cw = RD_RelativeDynamics(t, x_cw, F_cmd, T_cmd, params_cw);
    x_cw = x_cw + xdot_cw * dt_close;

    prop_remaining = prop_remaining - norm(F_cmd) * dt_close / (220*9.81);

    sim_log = append_log(sim_log, t, 3, 'Homing', x_cw(1:3), x_cw(4:6), ...
                     x_kep_chaser(1:3), x_kep_chaser(4:6), ...
                     [1;0;0;0], zeros(3,1), range, F_cmd, T_cmd, meas, prop_remaining);
end
t_now = t_now + (i-1)*dt_close;
fprintf('  Phase 3 complete. Range = %.1f m. t = %.0f s\n', range, t_now);

%% ----  PHASE 4: Closing  (200 m → 20 m) ------
fprintf('\n--- PHASE 4: Closing (200 m → 20 m) ---\n');
ph_cfg = cfg(5);

for i = 1:50000
    t = t_now + (i-1)*dt_close;
    range = norm(x_cw(1:3));
    if range < 20, break; end

    trueState = build_true_state(x_kep_chaser, x_kep_target, ...
                                  [1;0;0;0], zeros(3,1), zeros(3,1), ...
                                  sun_vec_ECI, mag_vec_ECI, ...
                                  0, false);
    trueState.rel_pos_LVLH = x_cw(1:3);
    trueState.rel_vel_LVLH = x_cw(4:6);

    meas = collect_measurements(ph_cfg, trueState);
    nav  = struct('rel_pos_LVLH', x_cw(1:3), 'rel_vel_LVLH', x_cw(4:6), ...
                  'range_m', range, 'att_quat', [1;0;0;0], 'ang_rate', zeros(3,1));

    [F_cmd, T_cmd] = RD_GuidanceControl(4, nav, struct(), ...
                      struct('n', n, 'mass_kg', mass_0, ...
                             'approach_axis', 'Vbar', ...
                             'approach_speed_ms', -0.05));

    params_cw = struct('mode','CW','n',n,'mass_kg',mass_0,'J2_enable',false);
    xdot_cw = RD_RelativeDynamics(t, x_cw, F_cmd, T_cmd, params_cw);
    x_cw = x_cw + xdot_cw * dt_close;

    prop_remaining = prop_remaining - norm(F_cmd) * dt_close / (220*9.81);

    sim_log = append_log(sim_log, t, 4, 'Closing', x_cw(1:3), x_cw(4:6), ...
                     x_kep_chaser(1:3), x_kep_chaser(4:6), ...
                     [1;0;0;0], zeros(3,1), range, F_cmd, T_cmd, meas, prop_remaining);
end
t_now = t_now + (i-1)*dt_close;
fprintf('  Phase 4 complete. Range = %.2f m. t = %.0f s\n', range, t_now);

%% ----  PHASE 5: Final Approach  (20 m → contact) ------
fprintf('\n--- PHASE 5: Final Approach (20 m → contact) ---\n');
ph_cfg = cfg(6);

for i = 1:50000
    t = t_now + (i-1)*dt_final;
    range = norm(x_cw(1:3));
    if range < 0.05, break; end  % ~5 cm – docking port contact

    trueState = build_true_state(x_kep_chaser, x_kep_target, ...
                                  [1;0;0;0], zeros(3,1), zeros(3,1), ...
                                  sun_vec_ECI, mag_vec_ECI, ...
                                  range, false);
    trueState.rel_pos_LVLH = x_cw(1:3);
    trueState.rel_vel_LVLH = x_cw(4:6);

    meas = collect_measurements(ph_cfg, trueState);
    nav  = struct('rel_pos_LVLH', x_cw(1:3), 'rel_vel_LVLH', x_cw(4:6), ...
                  'range_m', range, 'att_quat', [1;0;0;0], 'ang_rate', zeros(3,1));

    [F_cmd, T_cmd] = RD_GuidanceControl(5, nav, struct(), ...
                      struct('n', n, 'mass_kg', mass_0, ...
                             'approach_speed_ms', -0.01, ...
                             'q_dock_ref', [1;0;0;0]));

    params_cw = struct('mode','CW','n',n,'mass_kg',mass_0,'J2_enable',false);
    xdot_cw = RD_RelativeDynamics(t, x_cw, F_cmd, T_cmd, params_cw);
    x_cw = x_cw + xdot_cw * dt_final;

    prop_remaining = prop_remaining - norm(F_cmd) * dt_final / (220*9.81);

    sim_log = append_log(sim_log, t, 5, 'FinalApproach', x_cw(1:3), x_cw(4:6), ...
                     x_kep_chaser(1:3), x_kep_chaser(4:6), ...
                     [1;0;0;0], zeros(3,1), range, F_cmd, T_cmd, meas, prop_remaining);
end
t_now = t_now + (i-1)*dt_final;
fprintf('  Phase 5 complete. Range = %.4f m. t = %.0f s\n', range, t_now);

%% ----  PHASE 6: Docking  (contact + rigidisation) ------
fprintf('\n--- PHASE 6: Docking ---\n');
ph_cfg = cfg(7);

for i = 1:300  % 300 s max for docking
    t = t_now + (i-1)*dt_final;

    trueState = build_true_state(x_kep_chaser, x_kep_target, ...
                                  [1;0;0;0], zeros(3,1), zeros(3,1), ...
                                  sun_vec_ECI, mag_vec_ECI, ...
                                  0, true);  % contact_flag = true
    trueState.rel_pos_LVLH = x_cw(1:3);
    trueState.rel_vel_LVLH = x_cw(4:6);

    meas = collect_measurements(ph_cfg, trueState);
    nav  = struct('rel_pos_LVLH', x_cw(1:3), 'rel_vel_LVLH', x_cw(4:6), ...
                  'range_m', norm(x_cw(1:3)), 'att_quat', [1;0;0;0], ...
                  'ang_rate', zeros(3,1));

    [F_cmd, T_cmd] = RD_GuidanceControl(6, nav, struct(), ...
                      struct('n', n, 'mass_kg', mass_0));

    % During docking: spring-damper contact keeps spacecraft together
    x_cw(4:6) = x_cw(4:6) * 0.95;  % Damping of residual velocity

    sim_log = append_log(sim_log, t, 6, 'Docking', x_cw(1:3), x_cw(4:6), ...
                     x_kep_chaser(1:3), x_kep_chaser(4:6), ...
                     [1;0;0;0], zeros(3,1), norm(x_cw(1:3)), ...
                     F_cmd, T_cmd, meas, prop_remaining);
end
t_now = t_now + 300;
fprintf('  Phase 6 complete (Docked!). t = %.0f s  (%.1f hours)\n', ...
        t_now, t_now/3600);

%% -----------------------------------------------------------------------
%  5. Summary
% -----------------------------------------------------------------------
fprintf('\n==========================================================\n');
fprintf('  Mission Complete!\n');
fprintf('  Total simulated time : %.1f hours\n', t_now/3600);
fprintf('  Propellant used      : %.1f kg\n', 2000 - prop_remaining);
fprintf('  Total sim_log entries    : %d\n', numel(sim_log));
fprintf('==========================================================\n\n');

%% -----------------------------------------------------------------------
%  6. Post-processing plots
% -----------------------------------------------------------------------
RD_PlotResults(sim_log);

%% -----------------------------------------------------------------------
%  Local functions
% -----------------------------------------------------------------------

function sim_log = append_log(sim_log, t, ph_id, ph_name, rp, rv, pos, vel, q, w, r, F, T, meas, prop)
    entry.t              = t;
    entry.phase_id       = ph_id;
    entry.phase_name     = ph_name;
    entry.rel_pos_LVLH   = rp;
    entry.rel_vel_LVLH   = rv;
    entry.abs_pos_ECI    = pos;
    entry.abs_vel_ECI    = vel;
    entry.att_quat       = q;
    entry.ang_rate       = w;
    entry.range_m        = r;
    entry.F_cmd          = F;
    entry.T_cmd          = T;
    entry.sensor_meas    = meas;
    entry.propellant_kg  = prop;
    sim_log(end+1) = entry; %#ok<AGROW>
end

function meas = collect_measurements(ph_cfg, trueState)
%COLLECT_MEASUREMENTS  Collect measurements from all sensors in this phase
    meas = struct();
    for s = 1:numel(ph_cfg.sensors)
        stype = ph_cfg.sensors(s).type;
        try
            m = RD_SensorModels(stype, trueState, struct());
            meas.(stype) = m;
        catch
            meas.(stype) = struct('valid', false);
        end
    end
end

function ts = build_true_state(x_ch, x_tgt, q, w, rel_pos_hint, ...
                                sun_ECI, mag_ECI, gap_m, contact)
%BUILD_TRUE_STATE  Assemble the trueState struct for sensor models
    ts.pos_ECI      = x_ch(1:3);
    ts.vel_ECI      = x_ch(4:6);
    ts.att_quat     = q;
    ts.ang_rate     = w;
    ts.accel        = zeros(3,1);

    rel_pos_ECI     = x_ch(1:3) - x_tgt(1:3);
    rel_vel_ECI     = x_ch(4:6) - x_tgt(4:6);

    if norm(rel_pos_hint) > 0
        ts.rel_pos_LVLH = rel_pos_hint;
    else
        ts.rel_pos_LVLH = ECI_to_LVLH(rel_pos_ECI, x_tgt(1:3), x_tgt(4:6));
    end
    ts.rel_vel_LVLH = ECI_to_LVLH(rel_vel_ECI, x_tgt(1:3), x_tgt(4:6));
    ts.rel_att_quat = [1;0;0;0];
    ts.sun_vec_ECI  = sun_ECI;
    ts.mag_vec_body = mag_ECI;  % Simplified: ignore rotation
    ts.mag_vec_ECI  = mag_ECI;
    ts.gap_m        = gap_m;
    ts.contact_flag = contact;
end

function r_LVLH = ECI_to_LVLH(r_ECI, r_tgt, v_tgt)
%ECI_TO_LVLH  Convert ECI relative vector to LVLH frame
%   x = radial, y = along-track, z = cross-track
    r_hat = r_tgt / norm(r_tgt);
    h_vec = cross(r_tgt, v_tgt);
    z_hat = h_vec / norm(h_vec);
    y_hat = cross(z_hat, r_hat);
    R = [r_hat, y_hat, z_hat]';  % 3x3 rotation: ECI → LVLH
    r_LVLH = R * r_ECI;
end

function x_new = rk4(x, F, T, t, dt, params)
%RK4  Classic 4th-order Runge-Kutta integration
    k1 = RD_RelativeDynamics(t,       x,          F, T, params);
    k2 = RD_RelativeDynamics(t+dt/2,  x+dt/2*k1,  F, T, params);
    k3 = RD_RelativeDynamics(t+dt/2,  x+dt/2*k2,  F, T, params);
    k4 = RD_RelativeDynamics(t+dt,    x+dt*k3,    F, T, params);
    x_new = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end

function v = normalise(v)
    n = norm(v);
    if n > 1e-12, v = v/n; end
end

function print_phase_summary(cfg)
    fprintf('%-5s %-30s %-15s %-25s %-25s\n', ...
        'Phase', 'Name', 'Range [km]', 'Primary Sensors', 'Actuators');
    fprintf('%s\n', repmat('-',1,105));
    for p = 1:numel(cfg)
        ph = cfg(p);
        psen = {};
        for s = 1:numel(ph.sensors)
            if ph.sensors(s).primary, psen{end+1} = ph.sensors(s).type; end %#ok<AGROW>
        end
        acts = {ph.actuators.type};
        rng = sprintf('[%.0f-%.0f]', ph.range_km(1), ph.range_km(2));
        if isinf(ph.range_km(2)), rng = sprintf('>%.0f',ph.range_km(1)); end
        if ph.range_km(1)==0 && ph.range_km(2)==0, rng='Contact'; end
        fprintf('%-5d %-30s %-15s %-25s %-25s\n', ph.id, ph.name, rng, ...
            strjoin(psen, ','), strjoin(acts, ','));
    end
    fprintf('\n');
end
