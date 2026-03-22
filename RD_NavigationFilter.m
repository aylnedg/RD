function [x_est, P, nav] = RD_NavigationFilter(phase_id, x_est_prev, P_prev, ...
                                                  measurements, dt, params)
%RD_NAVIGATIONFILTER  Extended Kalman Filter (EKF) for R&D navigation
%   Fuses all available sensor measurements depending on the mission phase.
%
%   [x_est, P, nav] = RD_NavigationFilter(phase_id, x_est_prev, P_prev,
%                                          measurements, dt, params)
%
%   Navigation modes by phase (Fehse Ch.4):
%   Phase 0 – Launch  : INS + GPS absolute  (6-state: pos/vel ECI)
%   Phase 1 – Phasing : GPS + Rel.GPS + StarTracker  (12-state: abs+rel pos/vel)
%   Phase 2 – Far R&V : GPS + Radar + StarTracker    (12-state)
%   Phase 3 – Homing  : Radar + LIDAR + StarTracker  (12-state CW)
%   Phase 4 – Closing : LIDAR + VBS + StarTracker    (12-state CW)
%   Phase 5 – Final   : Flash LIDAR + VBS_HR + IMU   (12-state + attitude 4)
%   Phase 6 – Docking : DockingSensors + IMU          (6-state pose)
%
%   State vector (12-state relative nav, phases 1-4):
%     x = [rel_pos_LVLH(3); rel_vel_LVLH(3);
%           abs_pos_ECI(3);  abs_vel_ECI(3)]
%
%   State vector (16-state full nav, phases 5-6):
%     x = [rel_pos_LVLH(3); rel_vel_LVLH(3);
%           att_quat(4);     ang_rate(3);
%           gyro_bias(3)]
%
%   Inputs:
%     phase_id      integer  Current phase (0–6)
%     x_est_prev    [Nx1]    Previous state estimate
%     P_prev        [NxN]    Previous covariance matrix
%     measurements  struct   Struct with sensor measurement sub-structs
%                            (fields: .GPS, .RelGPS, .Radar, .LIDAR,
%                                     .VBS, .StarTracker, .IMU, etc.)
%     dt            scalar   Time step [s]
%     params        struct   Filter parameters (noise covariances, etc.)
%
%   Outputs:
%     x_est  [Nx1]   Updated state estimate
%     P      [NxN]   Updated covariance
%     nav    struct  Derived navigation quantities for GNC
%
% Reference:
%   Fehse, W. (2003). Automated Rendezvous and Docking of Spacecraft.
%   Cambridge University Press. Ch. 4.
%   Bar-Shalom, Y., Li, X.R., Kirubarajan, T. (2001). Estimation with
%   Applications to Tracking and Navigation. Wiley.
%
% See also: RD_SensorModels, RD_GuidanceControl, RD_RunAllPhases

% Copyright 2024 - Extended R&D Simulation

%% Determine state dimension from phase
if phase_id <= 2
    N = 12;   % [rel_pos(3); rel_vel(3); abs_pos(3); abs_vel(3)]
elseif phase_id <= 4
    N = 12;   % [rel_pos(3); rel_vel(3); abs_pos(3); abs_vel(3)]
elseif phase_id == 5
    N = 13;   % [rel_pos(3); rel_vel(3); att_quat(4); gyro_bias(3)]
else
    N = 6;    % [rel_pos_dp(3); rel_att_euler(3)]  (docking phase)
end

%% Initialise defaults if prev state size doesn't match
if isempty(x_est_prev) || numel(x_est_prev) ~= N
    x_est_prev = zeros(N, 1);
    if phase_id == 5
        x_est_prev(7) = 1;  % Identity quaternion scalar part
    end
end
if isempty(P_prev) || ~isequal(size(P_prev), [N N])
    P_prev = eye(N) * 1e6;
end

%% Mean motion for CW dynamics
n = getparam(params, 'n', 1.1e-3);   % rad/s

%% -----------------------------------------------------------------------
%  PREDICTION STEP (time update)
% -----------------------------------------------------------------------

if phase_id <= 4
    % Use CW state transition matrix for relative motion (phases 1-4)
    Phi_cw = CW_STM(n, dt);  % 6x6
    Phi    = blkdiag(Phi_cw, eye(6));  % 12x12 (rel + abs)

    % Process noise covariance Q
    q_rel = getparam(params, 'q_rel_pos', 1e-4);    % m²/s⁴
    q_abs = getparam(params, 'q_abs_pos', 1e-6);    % m²/s⁴
    Q_cw  = q_rel * [dt^3/3*eye(3), dt^2/2*eye(3); ...
                      dt^2/2*eye(3), dt*eye(3)];
    Q_abs = q_abs * [dt^3/3*eye(3), dt^2/2*eye(3); ...
                      dt^2/2*eye(3), dt*eye(3)];
    Q     = blkdiag(Q_cw, Q_abs);

    % Predicted state and covariance
    x_pred = Phi * x_est_prev;
    P_pred = Phi * P_prev * Phi' + Q;

elseif phase_id == 5
    % 6DOF navigation (relative pose + attitude + gyro bias)
    % State: [rel_pos(3), rel_vel(3), att_q(4), gyro_bias(3)]
    q0       = x_est_prev(7:10);
    omega_b  = x_est_prev(11:13);
    bias     = zeros(3,1);  % Could be estimated in x if augmented

    if isfield(measurements, 'IMU') && measurements.IMU.valid
        omega_meas = measurements.IMU.ang_rate;
    else
        omega_meas = omega_b;
    end

    % Attitude kinematics (quaternion integration, 1st order)
    q_dot = 0.5 * quatOmega(omega_meas) * q0;
    q_new = q0 + q_dot * dt;
    q_new = q_new / norm(q_new);

    % CW relative motion for position
    Phi_cw = CW_STM(n, dt);
    rp_new = Phi_cw * x_est_prev(1:6);

    x_pred = [rp_new; q_new; omega_meas];

    % Process noise
    q_rel  = getparam(params, 'q_rel_pos', 1e-4);
    q_att  = getparam(params, 'q_att',     1e-8);
    q_rate = getparam(params, 'q_rate',    1e-7);
    Q_6   = q_rel * [dt^3/3*eye(3), dt^2/2*eye(3); ...
                      dt^2/2*eye(3), dt*eye(3)];
    Q_att  = q_att  * eye(4) * dt;
    Q_rate = q_rate * eye(3) * dt;
    Q      = blkdiag(Q_6, Q_att, Q_rate);

    % Linearised F matrix for EKF (Jacobian of state transition)
    F_cw = CW_F(n);  % 6x6 continuous-time matrix
    F_att = -0.5 * skew(omega_meas);  % 3x3 attitude dynamics (simplified)
    F_lin = blkdiag(F_cw, F_att, zeros(3), zeros(3));
    Phi_lin = eye(N) + F_lin * dt;  % 1st order linearisation

    P_pred = Phi_lin * P_prev * Phi_lin' + Q;

else
    % Phase 6 docking: simplified 6-DOF pose estimate
    x_pred = x_est_prev;  % Near-static during docking
    P_pred = P_prev + eye(N) * getparam(params, 'q_dock', 1e-4) * dt;
end

%% -----------------------------------------------------------------------
%  UPDATE STEP (measurement update) – per sensor
% -----------------------------------------------------------------------

x_est = x_pred;
P     = P_pred;

%% GPS absolute position/velocity update (phases 0-3)
if phase_id <= 3 && isfield(measurements, 'GPS') && measurements.GPS.valid
    R_gps_p = getparam(params, 'R_GPS_pos_m2',  100) * eye(3);   % 10m sigma
    R_gps_v = getparam(params, 'R_GPS_vel_m2s', 0.01) * eye(3);  % 0.1 m/s sigma
    R_gps   = blkdiag(R_gps_p, R_gps_v);

    z_gps = [measurements.GPS.pos_ECI; measurements.GPS.vel_ECI];
    H_gps = [zeros(6,6), eye(6)];  % observes [abs_pos, abs_vel] part of x

    [x_est, P] = EKF_update(x_est, P, z_gps, H_gps, R_gps);
end

%% Relative GPS update (phases 1-3)
if phase_id >= 1 && phase_id <= 3 && ...
   isfield(measurements, 'RelativeGPS') && measurements.RelativeGPS.valid

    R_rgps_p = getparam(params, 'R_RGPS_pos_m2', 25)   * eye(3);  % 5m sigma
    R_rgps_v = getparam(params, 'R_RGPS_vel_m2s', 0.0025) * eye(3);
    R_rgps   = blkdiag(R_rgps_p, R_rgps_v);

    z_rgps = [measurements.RelativeGPS.rel_pos_LVLH; ...
              measurements.RelativeGPS.rel_vel_LVLH];
    H_rgps = [eye(6), zeros(6)];  % observes [rel_pos, rel_vel]

    [x_est, P] = EKF_update(x_est, P, z_rgps, H_rgps, R_rgps);
end

%% Radar update: range, range-rate, azimuth, elevation (phases 2,3)
if phase_id >= 2 && phase_id <= 3 && ...
   isfield(measurements, 'Radar') && measurements.Radar.valid

    R_radar = diag([measurements.Radar.valid * 400, ...  % 20m sigma²
                    0.04, 0.0025, 0.0025]);               % rr, az, el

    rel_pos_pred = x_est(1:3);
    r_pred    = norm(rel_pos_pred);
    rhat_pred = rel_pos_pred / max(r_pred, 1e-3);

    % Predicted measurements
    h_pred = [r_pred; ...
              dot(x_est(4:6), rhat_pred); ...
              atan2(rel_pos_pred(2), rel_pos_pred(1)); ...
              asin(rel_pos_pred(3) / max(r_pred, 1e-3))];

    % Measurement vector
    z_radar = [measurements.Radar.range_m; ...
               measurements.Radar.range_rate_ms; ...
               measurements.Radar.az_rad; ...
               measurements.Radar.el_rad];

    % Jacobian H (4x12, linearised around current estimate)
    H_r = jacobian_radar(rel_pos_pred, x_est(4:6));  % 4x6
    H_radar = [H_r, zeros(4,6)];  % append zeros for abs state part

    [x_est, P] = EKF_update_nonlinear(x_est, P, z_radar, h_pred, H_radar, R_radar);
end

%% LIDAR update: range, LOS (phases 3,4)
if phase_id >= 3 && phase_id <= 5 && ...
   isfield(measurements, 'LIDAR') && measurements.LIDAR.valid

    R_lidar = diag([1.0, 25e-6, 25e-6, 25e-6]);  % range:1m, LOS:5mrad sigma

    rel_pos_pred = x_est(1:3);
    r_pred    = norm(rel_pos_pred);
    rhat_pred = rel_pos_pred / max(r_pred, 1e-3);

    h_pred  = [r_pred; rhat_pred];
    z_lidar = [measurements.LIDAR.range_m; measurements.LIDAR.LOS_vec];

    H_r   = [rhat_pred', zeros(1,3); ...
             eye(3)/max(r_pred,1e-3) - rhat_pred*rhat_pred'/max(r_pred,1e-3), zeros(3)];
    if size(x_est,1) == 12
        H_lidar = [H_r, zeros(4,6)];
    else
        H_lidar = [H_r, zeros(4, size(x_est,1)-6)];
    end

    [x_est, P] = EKF_update_nonlinear(x_est, P, z_lidar, h_pred, H_lidar, R_lidar);
end

%% VBS Camera update: relative position (phases 3,4,5)
if phase_id >= 3 && phase_id <= 5 && ...
   isfield(measurements, 'VBS') && measurements.VBS.valid

    if phase_id == 5
        R_vbs = eye(3) * 0.0025;   % 0.05 m sigma (high-res)
    else
        R_vbs = eye(3) * 4.0;      % 2 m sigma (std camera)
    end

    z_vbs = measurements.VBS.rel_pos_m;
    H_vbs_rel = [eye(3), zeros(3)];  % observes rel_pos
    if size(x_est,1) == 12
        H_vbs = [H_vbs_rel, zeros(3,6)];
    else
        H_vbs = [H_vbs_rel, zeros(3, size(x_est,1)-6)];
    end

    [x_est, P] = EKF_update(x_est, P, z_vbs, H_vbs, R_vbs);
end

%% Star Tracker attitude update (phases 0-5)
if phase_id <= 5 && isfield(measurements, 'StarTracker') && ...
   measurements.StarTracker.valid && phase_id == 5

    % In 6DOF phase update attitude part of state
    R_st = (2.42e-5)^2 * eye(4);   % 5 arcsec crossbore sigma
    z_st = measurements.StarTracker.att_quat;
    H_st = [zeros(4,6), eye(4), zeros(4,3)];  % att_quat is states 7-10

    [x_est, P] = EKF_update(x_est, P, z_st, H_st, R_st);

    % Re-normalise quaternion
    x_est(7:10) = x_est(7:10) / norm(x_est(7:10));
end

%% Docking sensors update (phase 6)
if phase_id == 6 && isfield(measurements, 'DockingSensors') && ...
   measurements.DockingSensors.valid

    R_dock = diag([25e-6, 25e-6, 4e-6, 4e-6, 4e-6]);  % 5mm, 2mrad
    z_dock = [measurements.DockingSensors.alignment_m; ...
              measurements.DockingSensors.att_error_rad(2:4)];
    H_dock = eye(5,6);

    [x_est, P] = EKF_update(x_est, P, z_dock, H_dock, R_dock);
end

%% Symmetrise covariance
P = (P + P') / 2;

%% Build nav output struct for GNC
nav.rel_pos_LVLH = x_est(1:3);
nav.rel_vel_LVLH = x_est(4:6);
nav.range_m      = norm(x_est(1:3));
nav.range_rate_ms = dot(x_est(4:6), x_est(1:3)) / max(norm(x_est(1:3)), 1e-3);
nav.P_rel_pos    = P(1:3, 1:3);
nav.sigma_range  = sqrt(nav.P_rel_pos(1,1) + nav.P_rel_pos(2,2) + nav.P_rel_pos(3,3));

if phase_id == 5 && N == 13
    nav.att_quat  = x_est(7:10) / norm(x_est(7:10));
    nav.ang_rate  = x_est(11:13);
elseif phase_id <= 4 && isfield(measurements, 'StarTracker') && ...
       measurements.StarTracker.valid
    nav.att_quat  = measurements.StarTracker.att_quat;
    nav.ang_rate  = zeros(3,1);
else
    nav.att_quat  = [1;0;0;0];
    nav.ang_rate  = zeros(3,1);
end

end

%% -----------------------------------------------------------------------
%  EKF Update Helpers
% -----------------------------------------------------------------------

function [x_upd, P_upd] = EKF_update(x_pred, P_pred, z, H, R)
%EKF_UPDATE  Linear EKF measurement update
    S = H * P_pred * H' + R;
    K = P_pred * H' / S;
    innov = z - H * x_pred;
    x_upd = x_pred + K * innov;
    n = size(P_pred, 1);
    P_upd = (eye(n) - K * H) * P_pred;
end

function [x_upd, P_upd] = EKF_update_nonlinear(x_pred, P_pred, z, h, H, R)
%EKF_UPDATE_NONLINEAR  EKF update using nonlinear measurement (innovation = z - h)
    S = H * P_pred * H' + R;
    K = P_pred * H' / S;
    innov = z - h;
    % Angle wrapping for az/el
    innov(3:end) = angwrap(innov(3:end));
    x_upd = x_pred + K * innov;
    n = size(P_pred, 1);
    P_upd = (eye(n) - K * H) * P_pred;
end

%% -----------------------------------------------------------------------
%  CW Matrix Helpers
% -----------------------------------------------------------------------

function Phi = CW_STM(n, dt)
%CW_STM  6x6 CW state transition matrix  (Fehse Eq. 3.42)
    s = sin(n*dt); c = cos(n*dt);
    Phi = [4-3*c,      0, 0,  s/n,        2*(1-c)/n,       0;   ...
           6*(s-n*dt), 1, 0, -2*(1-c)/n, (4*s-3*n*dt)/n,  0;   ...
           0,          0, c,  0,          0,               s/n; ...
           3*n*s,      0, 0,  c,          2*s,             0;   ...
          -6*n*(1-c),  0, 0, -2*s,        4*c-3,           0;   ...
           0,          0,-n*s, 0,         0,               c];
end

function F = CW_F(n)
%CW_F  Continuous-time CW dynamics matrix (6x6)
    F = [0, 0, 0, 1, 0, 0;   ...
         0, 0, 0, 0, 1, 0;   ...
         0, 0, 0, 0, 0, 1;   ...
         3*n^2, 0, 0, 0, 2*n, 0;  ...
         0,     0, 0,-2*n, 0, 0;  ...
         0,     0,-n^2, 0, 0, 0];
end

function H = jacobian_radar(r, v)
%JACOBIAN_RADAR  4x6 Jacobian for [range; rangerate; az; el] w.r.t. [rel_pos; rel_vel]
    rn   = norm(r);
    rhat = r / max(rn, 1e-3);

    % drange/dr  drange/dv
    dR_dr = rhat';          dR_dv = zeros(1,3);
    % drangerate/dr  drangerate/dv
    dRd_dr = (v - dot(v,rhat)*rhat)' / max(rn, 1e-3);
    dRd_dv = rhat';
    % daz/dr (atan2 partial)
    daz_dr = [-r(2), r(1), 0] / max(r(1)^2 + r(2)^2, 1e-6);
    % del/dr (asin partial)
    dz_rn  = (rn^2 - r(3)^2);
    del_dr = [-r(1)*r(3), -r(2)*r(3), rn^2-r(3)^2] / max(rn^3 * sqrt(dz_rn), 1e-6);

    H = [dR_dr,  dR_dv;  ...
         dRd_dr, dRd_dv; ...
         daz_dr, zeros(1,3); ...
         del_dr, zeros(1,3)];
end

function Omega = quatOmega(w)
%QUATOPMEGA  4x4 quaternion rate matrix for q_dot = 0.5*Omega*q
    Omega = [ 0,    -w(1), -w(2), -w(3); ...
              w(1),  0,     w(3), -w(2); ...
              w(2), -w(3),  0,     w(1); ...
              w(3),  w(2), -w(1),  0  ];
end

function S = skew(v)
    S = [ 0,    -v(3),  v(2); ...
          v(3),  0,    -v(1); ...
         -v(2),  v(1),  0   ];
end

function a = angwrap(a)
    a = mod(a + pi, 2*pi) - pi;
end

function val = getparam(params, field, default)
    if isfield(params, field)
        val = params.(field);
    else
        val = default;
    end
end
