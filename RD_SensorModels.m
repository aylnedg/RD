function meas = RD_SensorModels(sensorType, trueState, params)
%RD_SENSORMODELS  Sensor measurement models for all R&D phases
%   meas = RD_SensorModels(sensorType, trueState, params)
%   returns a noisy measurement struct for the requested sensor type.
%
%   SENSOR TYPES (sensorType string):
%     'GPS'                  - Absolute position/velocity (ECI)
%     'RelativeGPS'          - Relative position/velocity (LVLH)
%     'Radar_RF'             - Range, range-rate, azimuth, elevation
%     'LIDAR_Scanning'       - Range, bearing, LOS vector
%     'LIDAR_Flash'          - Range, range-rate, 3-D point cloud
%     'VBS_Camera'           - Relative pose from feature tracking
%     'VBS_Camera_HighRes'   - High-res pose + docking port
%     'StarTracker'          - Attitude quaternion (ICRF)
%     'SunSensor_Coarse'     - Coarse sun direction (body frame)
%     'SunSensor_Fine'       - Fine sun direction (body frame)
%     'EarthHorizonSensor'   - Nadir direction (body frame)
%     'IMU'                  - Angular rate + specific force
%     'Magnetometer'         - Magnetic field vector (body frame)
%     'ProximitySensors'     - Gap distance + contact flag
%     'DockingSensors'       - Alignment + capture flag
%     'ContactForceSensors'  - Force/torque at contact
%
%   trueState : struct with fields relevant to the sensor:
%     .pos_ECI      [3x1] m    chaser position in ECI
%     .vel_ECI      [3x1] m/s  chaser velocity in ECI
%     .att_quat     [4x1]      scalar-first body->ICRF quaternion
%     .ang_rate     [3x1] rad/s body angular rate
%     .accel        [3x1] m/s² specific force (body frame)
%     .rel_pos_LVLH [3x1] m    relative position in LVLH
%     .rel_vel_LVLH [3x1] m/s  relative velocity in LVLH
%     .rel_att_quat [4x1]      relative attitude quat
%     .sun_vec_ECI  [3x1]      unit sun direction in ECI
%     .mag_vec_ECI  [3x1] T    Earth magnetic field in ECI
%     .gap_m        scalar m   docking port gap
%     .contact_flag logical    docking contact flag
%
%   params : (optional) struct to override default noise parameters.
%
% Reference:
%   Fehse, W. (2003). Automated Rendezvous and Docking of Spacecraft.
%   Cambridge University Press. ISBN 0-521-82492-3.
%
% See also: RD_Fehse_PhaseConfig, RD_NavigationFilter

% Copyright 2024 - Extended R&D Simulation

if nargin < 3 || isempty(params)
    params = struct();
end

meas = struct();

switch sensorType

    % ===================================================================
    case 'GPS'
    % -------------------------------------------------------------------
    % Absolute GPS receiver: position and velocity in ECI frame
    % Typical: SEM-E GPS, accuracy ~10 m (pos) / 0.1 m/s (vel) 1-sigma
    % Active in: Phases 0,1,2,3
    % -------------------------------------------------------------------
        sigma_pos = getparam(params, 'sigma_pos_m',   10);   % m
        sigma_vel = getparam(params, 'sigma_vel_ms',  0.1);  % m/s

        meas.pos_ECI = trueState.pos_ECI + sigma_pos * randn(3,1);
        meas.vel_ECI = trueState.vel_ECI + sigma_vel * randn(3,1);
        meas.valid   = true;

    % ===================================================================
    case 'RelativeGPS'
    % -------------------------------------------------------------------
    % Differential / carrier-phase GPS for relative navigation
    % Accuracy: ~5 m (pos), 0.05 m/s (vel) 1-sigma for baselines < 100 km
    % Active in: Phases 1,2,3
    % -------------------------------------------------------------------
        sigma_rpos = getparam(params, 'sigma_rel_pos_m',  5);    % m
        sigma_rvel = getparam(params, 'sigma_rel_vel_ms', 0.05); % m/s

        meas.rel_pos_LVLH = trueState.rel_pos_LVLH + sigma_rpos * randn(3,1);
        meas.rel_vel_LVLH = trueState.rel_vel_LVLH + sigma_rvel * randn(3,1);
        meas.valid = true;

    % ===================================================================
    case 'Radar_RF'
    % -------------------------------------------------------------------
    % Rendezvous radar (RF): range, range-rate, azimuth, elevation
    % Typical: VHF/UHF Kurs or S-band radar
    % Accuracy: range ~20 m, range-rate ~0.2 m/s, angles ~50 mrad
    % Active in: Phases 2,3
    % -------------------------------------------------------------------
        sigma_r    = getparam(params, 'sigma_range_m',      20);   % m
        sigma_rdot = getparam(params, 'sigma_rangerate_ms', 0.2);  % m/s
        sigma_ang  = getparam(params, 'sigma_angle_rad',    0.05); % rad

        relvec  = trueState.rel_pos_LVLH;  % LVLH, x=radial, y=along-track, z=cross-track
        range   = norm(relvec);
        relvec_u = relvec / max(range, 1e-3);

        range_rate = dot(trueState.rel_vel_LVLH, relvec_u);
        az_true    = atan2(relvec(2), relvec(1));
        el_true    = asin(relvec(3) / max(range, 1e-3));

        meas.range_m       = range      + sigma_r    * randn();
        meas.range_rate_ms = range_rate + sigma_rdot * randn();
        meas.az_rad        = az_true    + sigma_ang  * randn();
        meas.el_rad        = el_true    + sigma_ang  * randn();
        meas.valid = (range < 50000);  % 50 km max range

    % ===================================================================
    case 'LIDAR_Scanning'
    % -------------------------------------------------------------------
    % Scanning LIDAR (e.g., Rendezvous Laser Radar RLR):
    % Provides range, LOS direction, bearing
    % Accuracy: range ~1 m, bearing ~10 mrad, LOS ~5 mrad
    % Active in: Phases 3,4
    % -------------------------------------------------------------------
        sigma_r   = getparam(params, 'sigma_range_m',   1.0);  % m
        sigma_los = getparam(params, 'sigma_LOS_rad',   0.005);% rad

        relvec = trueState.rel_pos_LVLH;
        range  = norm(relvec);
        los_true = relvec / max(range, 1e-3);

        meas.range_m   = range + sigma_r   * randn();
        meas.LOS_vec   = normalise(los_true + sigma_los * randn(3,1));
        meas.bearing_rad = acos(min(max(dot([1;0;0], los_true), -1), 1));
        meas.valid = (range < 5000 && range > 0.1);

    % ===================================================================
    case 'LIDAR_Flash'
    % -------------------------------------------------------------------
    % Flash LIDAR / 3-D imaging LIDAR:
    % Provides full 3-D point cloud, range, range-rate, relative attitude
    % Accuracy: range ~0.2 m, range-rate ~0.02 m/s, att ~0.05 deg
    % Active in: Phases 4,5
    % -------------------------------------------------------------------
        sigma_r    = getparam(params, 'sigma_range_m',     0.2);   % m
        sigma_rdot = getparam(params, 'sigma_rangerate_ms',0.02);  % m/s
        sigma_att  = getparam(params, 'sigma_att_rad',     0.001); % rad (~0.05 deg)

        relvec  = trueState.rel_pos_LVLH;
        range   = norm(relvec);
        relvec_u = relvec / max(range, 1e-3);
        range_rate = dot(trueState.rel_vel_LVLH, relvec_u);

        meas.range_m       = range      + sigma_r    * randn();
        meas.range_rate_ms = range_rate + sigma_rdot * randn();
        meas.rel_att_quat  = addAttNoise(trueState.rel_att_quat, sigma_att);
        meas.point_cloud   = generate_mock_point_cloud(trueState, 64);
        meas.valid = (range < 300 && range > 0.01);

    % ===================================================================
    case 'VBS_Camera'
    % -------------------------------------------------------------------
    % Vision-Based Sensor (monocular/stereo camera):
    % Relative position, attitude, LOS from feature matching
    % Accuracy: rel_pos ~2 m, att ~0.5 deg, LOS ~3 mrad
    % Active in: Phases 3,4
    % -------------------------------------------------------------------
        sigma_rpos = getparam(params, 'sigma_rel_pos_m', 2.0);  % m
        sigma_att  = getparam(params, 'sigma_att_rad',   0.009);% rad (~0.5 deg)
        sigma_los  = getparam(params, 'sigma_LOS_rad',   0.003);% rad

        relvec = trueState.rel_pos_LVLH;
        range  = norm(relvec);
        los_true = relvec / max(range, 1e-3);

        meas.rel_pos_m    = trueState.rel_pos_LVLH + sigma_rpos * randn(3,1);
        meas.rel_att_quat = addAttNoise(trueState.rel_att_quat, sigma_att);
        meas.LOS_vec      = normalise(los_true + sigma_los * randn(3,1));
        meas.valid = (range < 5000 && range > 1.0);

    % ===================================================================
    case 'VBS_Camera_HighRes'
    % -------------------------------------------------------------------
    % High-resolution VBS for final approach / docking:
    % Docking port pose with sub-cm position, sub-mrad attitude
    % Accuracy: rel_pos ~0.05 m, att ~0.02 deg, LOS ~0.5 mrad
    % Active in: Phases 5,6
    % -------------------------------------------------------------------
        sigma_rpos = getparam(params, 'sigma_rel_pos_m', 0.05);   % m
        sigma_att  = getparam(params, 'sigma_att_rad',   0.00035);% rad (~0.02 deg)
        sigma_los  = getparam(params, 'sigma_LOS_rad',   0.0005); % rad

        relvec = trueState.rel_pos_LVLH;
        range  = norm(relvec);
        los_true = relvec / max(range, 1e-3);

        meas.rel_pos_m        = trueState.rel_pos_LVLH + sigma_rpos * randn(3,1);
        meas.rel_att_quat     = addAttNoise(trueState.rel_att_quat, sigma_att);
        meas.LOS_vec          = normalise(los_true + sigma_los * randn(3,1));
        meas.docking_port_pos = trueState.rel_pos_LVLH + sigma_rpos * randn(3,1);
        meas.valid = (range < 100);

    % ===================================================================
    case 'StarTracker'
    % -------------------------------------------------------------------
    % Star tracker: absolute inertial attitude quaternion
    % Typical: Terma HE-5AS, accuracy ~5 arcsec cross-bore, ~30 arcsec roll
    % Active in: All phases
    % -------------------------------------------------------------------
        sigma_cb  = getparam(params, 'sigma_crossbore_rad', 2.42e-5); % 5 arcsec
        sigma_rol = getparam(params, 'sigma_roll_rad',      1.45e-4); % 30 arcsec

        % Cross-bore errors affect x,y axes; roll error affects z-axis
        err_vec = [sigma_cb * randn(2,1); sigma_rol * randn()];
        meas.att_quat = addAttNoise(trueState.att_quat, norm(err_vec));
        meas.valid = true;  % Simplified - ignore sun/moon exclusion angles

    % ===================================================================
    case 'SunSensor_Coarse'
    % -------------------------------------------------------------------
    % Coarse sun sensor (CSS): 2-axis sun direction in body frame
    % Accuracy: ~1 deg 1-sigma
    % Active in: Phases 0,1 (always available if not in eclipse)
    % -------------------------------------------------------------------
        sigma_ang = getparam(params, 'sigma_angle_rad', 0.01745); % 1 deg

        DCM = quat2dcm_sf(trueState.att_quat);  % ECI -> body
        sun_body = DCM * trueState.sun_vec_ECI;
        meas.sun_vec_body = normalise(sun_body + sigma_ang * randn(3,1));
        meas.valid = (dot(sun_body, sun_body) > 0.5);  % In sunlight

    % ===================================================================
    case 'SunSensor_Fine'
    % -------------------------------------------------------------------
    % Fine sun sensor (FSS / digital sun sensor):
    % Accuracy: ~0.05 deg 1-sigma
    % Active in: Phases 1,2,3
    % -------------------------------------------------------------------
        sigma_ang = getparam(params, 'sigma_angle_rad', 8.73e-4); % 0.05 deg

        DCM = quat2dcm_sf(trueState.att_quat);
        sun_body = DCM * trueState.sun_vec_ECI;
        meas.sun_vec_body = normalise(sun_body + sigma_ang * randn(3,1));
        meas.valid = (dot(sun_body, sun_body) > 0.5);

    % ===================================================================
    case 'EarthHorizonSensor'
    % -------------------------------------------------------------------
    % Earth horizon / infra-red earth sensor (IRES):
    % Provides nadir direction in body frame
    % Accuracy: ~0.1 deg pitch/roll 1-sigma
    % Active in: Phases 0,1,2
    % -------------------------------------------------------------------
        sigma_ang = getparam(params, 'sigma_angle_rad', 1.745e-3); % 0.1 deg

        DCM = quat2dcm_sf(trueState.att_quat);
        nadir_ECI = -trueState.pos_ECI / norm(trueState.pos_ECI);
        nadir_body = DCM * nadir_ECI;
        meas.nadir_vec = normalise(nadir_body + sigma_ang * randn(3,1));
        meas.valid = true;

    % ===================================================================
    case 'IMU'
    % -------------------------------------------------------------------
    % Inertial Measurement Unit: rate gyros + accelerometers
    % Typical: Honeywell GG1320 (ring-laser gyro)
    % Gyro: bias drift 0.01 deg/h, ARW 0.005 deg/sqrt(h)
    % Accel: bias 50 µg, noise 50 µg/sqrt(Hz)
    % Active in: All phases
    % -------------------------------------------------------------------
        % Gyroscope noise parameters
        arw       = getparam(params, 'ARW_degph',    0.005);  % deg/sqrt(h)
        bias_gyro = getparam(params, 'bias_dph',     0.01);   % deg/h
        dt        = getparam(params, 'dt_s',         0.0025); % 400 Hz

        sigma_gyro = (arw*pi/180/60) / sqrt(dt);  % convert to rad/s/sqrt(sample)
        bias_g_rads = (bias_gyro * pi/180 / 3600); % rad/s

        % Accelerometer noise parameters
        noise_accel = getparam(params, 'noise_ug',   50);     % µg/sqrt(Hz)
        bias_accel  = getparam(params, 'bias_ug',    50);     % µg

        sigma_accel = (noise_accel * 1e-6 * 9.81) / sqrt(dt);
        bias_a_ms2  = bias_accel * 1e-6 * 9.81;

        meas.ang_rate = trueState.ang_rate + bias_g_rads + sigma_gyro * randn(3,1);
        meas.accel    = trueState.accel    + bias_a_ms2  + sigma_accel * randn(3,1);
        meas.valid    = true;

    % ===================================================================
    case 'Magnetometer'
    % -------------------------------------------------------------------
    % Fluxgate magnetometer: Earth magnetic field in body frame
    % Accuracy: ~100 nT 1-sigma  (IGRF model accuracy ~200 nT)
    % Active in: Phases 0,1,2 (LEO only)
    % -------------------------------------------------------------------
        sigma_B = getparam(params, 'sigma_field_nT', 100); % nT

        DCM = quat2dcm_sf(trueState.att_quat);
        B_body = DCM * trueState.mag_vec_ECI;
        meas.mag_field_body = B_body + (sigma_B * 1e-9) * randn(3,1);
        meas.valid = true;

    % ===================================================================
    case 'ProximitySensors'
    % -------------------------------------------------------------------
    % Contact proximity sensors (micro-switch / capacitive):
    % Gap distance and contact flag for pre-contact detection
    % Accuracy: ~10 mm gap, contact flag (boolean)
    % Active in: Phase 5 (final approach close-in)
    % -------------------------------------------------------------------
        sigma_gap = getparam(params, 'sigma_gap_m', 0.01); % 10 mm

        meas.gap_m       = trueState.gap_m + sigma_gap * randn();
        meas.gap_m       = max(meas.gap_m, 0);
        meas.contact_flag = (trueState.gap_m <= 0);
        meas.valid = (trueState.gap_m < 0.5);  % 50 cm activation range

    % ===================================================================
    case 'DockingSensors'
    % -------------------------------------------------------------------
    % Docking mechanism sensors (alignment + capture ring sensors):
    % Misalignment and lateral/angular offset at docking interface
    % Accuracy: alignment ~5 mm, attitude ~2 mrad
    % Active in: Phase 6 (docking)
    % -------------------------------------------------------------------
        sigma_align = getparam(params, 'sigma_align_m',  0.005); % 5 mm
        sigma_att   = getparam(params, 'sigma_att_rad',  0.002); % 2 mrad

        meas.alignment_m    = trueState.rel_pos_LVLH(1:2) + sigma_align * randn(2,1);
        meas.att_error_rad  = addAttNoise(trueState.rel_att_quat, sigma_att);
        meas.capture_flag   = (norm(trueState.rel_pos_LVLH) < 0.1) && ...
                              trueState.contact_flag;
        meas.valid = true;

    % ===================================================================
    case 'ContactForceSensors'
    % -------------------------------------------------------------------
    % 6-DOF force/torque sensor at docking interface:
    % Measures contact forces during capture and rigidisation
    % Accuracy: force ~0.5 N, torque ~0.1 Nm
    % Active in: Phase 6 (docking, after contact)
    % -------------------------------------------------------------------
        sigma_F  = getparam(params, 'sigma_force_N',  0.5); % N
        sigma_T  = getparam(params, 'sigma_torque_Nm',0.1); % Nm

        if isfield(trueState, 'contact_force_N')
            meas.force_N   = trueState.contact_force_N  + sigma_F * randn(3,1);
            meas.torque_Nm = trueState.contact_torque_Nm+ sigma_T * randn(3,1);
        else
            meas.force_N   = sigma_F * randn(3,1);
            meas.torque_Nm = sigma_T * randn(3,1);
        end
        meas.valid = trueState.contact_flag;

    otherwise
        error('RD_SensorModels: Unknown sensor type ''%s''', sensorType);
end

end

%% -----------------------------------------------------------------------
%  Local helper functions
% -----------------------------------------------------------------------

function val = getparam(params, field, default)
%GETPARAM  Return params.field if it exists, otherwise return default
    if isfield(params, field)
        val = params.(field);
    else
        val = default;
    end
end

function q_noisy = addAttNoise(q_true, sigma_rad)
%ADDATTNOISE  Add angular noise to a scalar-first quaternion
%   Generates a small rotation of magnitude sigma_rad about a random axis
    axis = randn(3,1);
    if norm(axis) < 1e-10
        axis = [1;0;0];
    end
    axis = axis / norm(axis);
    angle = sigma_rad * randn();
    dq = [cos(angle/2); sin(angle/2) * axis];  % scalar-first perturbation
    % Quaternion product (scalar-first convention): q_noisy = dq ⊗ q_true
    q_noisy = quatmultiply_sf(dq, q_true);
    q_noisy = q_noisy / norm(q_noisy);
end

function q_ab = quatmultiply_sf(q_a, q_b)
%QUATMULTIPLY_SF  Quaternion product, scalar-first convention
    w1 = q_a(1); v1 = q_a(2:4);
    w2 = q_b(1); v2 = q_b(2:4);
    q_ab = [w1*w2 - dot(v1,v2); ...
            w1*v2 + w2*v1 + cross(v1,v2)];
end

function DCM = quat2dcm_sf(q)
%QUAT2DCM_SF  Rotation matrix from scalar-first quaternion (body <- ECI)
    w = q(1); x = q(2); y = q(3); z = q(4);
    DCM = [1-2*(y^2+z^2),   2*(x*y+w*z),   2*(x*z-w*y); ...
           2*(x*y-w*z),   1-2*(x^2+z^2),   2*(y*z+w*x); ...
           2*(x*z+w*y),   2*(y*z-w*x),   1-2*(x^2+y^2)];
end

function v_u = normalise(v)
%NORMALISE  Return unit vector, guarded against zero-norm
    n = norm(v);
    if n < 1e-12
        v_u = v;
    else
        v_u = v / n;
    end
end

function cloud = generate_mock_point_cloud(trueState, N)
%GENERATE_MOCK_POINT_CLOUD  Returns N 3-D points representing target surface
%   Simplified: uniform scatter on sphere of 5 m radius around target
    r_target = trueState.rel_pos_LVLH;
    pts = randn(3, N);
    pts = pts ./ vecnorm(pts);       % unit vectors
    pts = pts * 5;                   % 5 m target radius
    cloud = pts + r_target;          % translate to target location
end
