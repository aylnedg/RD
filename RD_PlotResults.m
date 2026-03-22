function RD_PlotResults(sim_log)
%RD_PLOTRESULTS  Visualise the complete R&D mission simulation results
%
%   RD_PlotResults(sim_log)  generates a comprehensive set of plots showing:
%     1. Range vs Time (all phases, with phase boundaries)
%     2. Relative trajectory in LVLH frame (V-bar approach diagram)
%     3. Per-phase relative trajectory (CW plane + cross-track)
%     4. Commanded forces and propellant usage
%     5. Sensor availability matrix per phase
%
%   Input:
%     sim_log  struct array  Simulation sim_log from RD_RunAllPhases
%
% See also: RD_RunAllPhases, RD_Fehse_PhaseConfig

% Copyright 2024 - Extended R&D Simulation

if isempty(sim_log)
    warning('RD_PlotResults: empty sim_log'); return;
end

% Extract time-series
t_s      = [sim_log.t];
phase_id = [sim_log.phase_id];
range_m  = [sim_log.range_m];
rel_pos  = cell2mat({sim_log.rel_pos_LVLH});   % 3 × N
F_cmd    = cell2mat({sim_log.F_cmd});           % 3 × N
prop     = [sim_log.propellant_kg];

phase_colors = [0.8 0.2 0.2;   % 0 Launch – red
                0.2 0.6 0.9;   % 1 Phasing – blue
                0.2 0.8 0.4;   % 2 Far R&V – green
                0.9 0.7 0.1;   % 3 Homing – orange
                0.7 0.2 0.8;   % 4 Closing – purple
                0.2 0.9 0.9;   % 5 Final – cyan
                0.9 0.3 0.6];  % 6 Docking – pink

phase_names = {'Launch','Phasing','Far R&V','Homing','Closing', ...
               'Final Appr.','Docking'};

%% -----------------------------------------------------------------------
%  Figure 1: Range vs Time
% -----------------------------------------------------------------------
figure('Name','R&D Mission: Range vs Time','NumberTitle','off', ...
       'Position',[100 100 1200 500]);

semilogy(t_s/3600, range_m, 'k-', 'LineWidth', 1.5);
hold on; grid on;

% Shade phases
ax = gca; ylim_cur = get(ax,'YLim');
for ph = 0:6
    idx = find(phase_id == ph);
    if isempty(idx), continue; end
    t_start = t_s(idx(1))/3600;
    t_end   = t_s(idx(end))/3600;
    patch([t_start t_end t_end t_start], [1e-1 1e-1 1e6 1e6], ...
          phase_colors(ph+1,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    text((t_start+t_end)/2, max(range_m)*0.5, phase_names{ph+1}, ...
         'HorizontalAlignment','center','FontSize',8,'Rotation',90, ...
         'Color', phase_colors(ph+1,:)*0.7);
end

% Phase boundary markers
[~, trans_idx] = unique(phase_id, 'first');
for k = 1:numel(trans_idx)
    xline(t_s(trans_idx(k))/3600, '--', 'Color',[0.5 0.5 0.5], 'LineWidth',0.8);
end

xlabel('Mission Elapsed Time [hours]');
ylabel('Range to Target [m]');
title('Spacecraft R&D: Range vs Time – All Phases (Fehse 2003)');
ylim([0.01 max(range_m)*2]);

%% -----------------------------------------------------------------------
%  Figure 2: LVLH Relative Trajectory (V-bar approach view)
% -----------------------------------------------------------------------
figure('Name','R&D Mission: Relative Trajectory (LVLH)','NumberTitle','off', ...
       'Position',[100 200 1200 800]);

subplot(2,2,1);  % Full mission view
for ph = 0:6
    idx = find(phase_id == ph);
    if isempty(idx), continue; end
    x_km = rel_pos(2, idx)/1e3;  % Along-track (y)
    y_km = rel_pos(1, idx)/1e3;  % Radial (x)
    plot(x_km, y_km, '-', 'Color', phase_colors(ph+1,:), 'LineWidth', 1.5, ...
         'DisplayName', phase_names{ph+1});
    hold on;
end
plot(0, 0, 'k^', 'MarkerSize', 10, 'MarkerFaceColor','k', ...
     'DisplayName', 'Target (ISS)');
grid on; xlabel('Along-Track [km]'); ylabel('Radial [km]');
title('Full Mission Trajectory in LVLH'); legend('Location','best','FontSize',7);

subplot(2,2,2);  % Close-range view (phases 3-5)
for ph = 3:5
    idx = find(phase_id == ph);
    if isempty(idx), continue; end
    x_m = rel_pos(2, idx);   % Along-track [m]
    y_m = rel_pos(1, idx);   % Radial [m]
    plot(x_m, y_m, '-', 'Color', phase_colors(ph+1,:), 'LineWidth', 2, ...
         'DisplayName', phase_names{ph+1});
    hold on;
end
plot(0, 0, 'k^', 'MarkerSize', 10, 'MarkerFaceColor','k', ...
     'DisplayName', 'Target');
grid on; xlabel('Along-Track (V-bar) [m]'); ylabel('Radial (R-bar) [m]');
title('Close-Range: Phases 3-5 (Homing → Final Approach)');
legend('Location','best','FontSize',8);

subplot(2,2,3);  % Cross-track vs along-track
for ph = 0:6
    idx = find(phase_id == ph);
    if isempty(idx), continue; end
    x_m = rel_pos(2, idx);   % Along-track
    z_m = rel_pos(3, idx);   % Cross-track
    plot(x_m/1e3, z_m, '-', 'Color', phase_colors(ph+1,:), 'LineWidth', 1.5);
    hold on;
end
grid on; xlabel('Along-Track [km]'); ylabel('Cross-Track [m]');
title('Out-of-Plane Motion');

subplot(2,2,4);  % Approach corridor check (phases 4-5)
for ph = 4:5
    idx = find(phase_id == ph);
    if isempty(idx), continue; end
    rng = range_m(idx);
    lat = sqrt(rel_pos(1,idx).^2 + rel_pos(3,idx).^2);
    corridor = tand(15) * abs(rel_pos(2, idx));  % ±15° cone
    plot(rng, lat, '-', 'Color', phase_colors(ph+1,:), 'LineWidth', 1.5, ...
         'DisplayName', phase_names{ph+1});
    hold on;
    plot(rng, corridor, ':', 'Color', phase_colors(ph+1,:)*0.7, 'LineWidth', 1);
end
set(gca, 'XDir','reverse');
grid on; xlabel('Range [m]'); ylabel('Lateral Deviation [m]');
title('Approach Corridor Check (±15° cone)');
legend({'Phase 4 traj','Phase 4 corridor','Phase 5 traj','Phase 5 corridor'}, ...
       'Location','best','FontSize',7);

%% -----------------------------------------------------------------------
%  Figure 3: Forces and Propellant
% -----------------------------------------------------------------------
figure('Name','R&D Mission: Control Effort','NumberTitle','off', ...
       'Position',[150 300 1000 500]);

subplot(2,1,1);
F_mag = vecnorm(F_cmd);
for ph = 0:6
    idx = find(phase_id == ph);
    if isempty(idx), continue; end
    plot(t_s(idx)/3600, F_mag(idx), '-', 'Color', phase_colors(ph+1,:), ...
         'DisplayName', phase_names{ph+1});
    hold on;
end
grid on; xlabel('Time [hours]'); ylabel('|F_{cmd}| [N]');
title('Commanded Force Magnitude'); legend('Location','best','FontSize',7);

subplot(2,1,2);
plot(t_s/3600, prop, 'b-', 'LineWidth', 1.5);
grid on; xlabel('Time [hours]'); ylabel('Propellant Remaining [kg]');
title('Propellant Mass');

%% -----------------------------------------------------------------------
%  Figure 4: Sensor Availability Matrix
% -----------------------------------------------------------------------
cfg = RD_Fehse_PhaseConfig();

all_sensors = {};
for p = 1:numel(cfg)
    for s = 1:numel(cfg(p).sensors)
        stype = cfg(p).sensors(s).type;
        if ~any(strcmp(all_sensors, stype))
            all_sensors{end+1} = stype; %#ok<AGROW>
        end
    end
end

avail_matrix = zeros(numel(all_sensors), 7);  % sensors x phases
primary_matrix = zeros(numel(all_sensors), 7);

for p = 1:numel(cfg)
    for s = 1:numel(cfg(p).sensors)
        stype = cfg(p).sensors(s).type;
        row = find(strcmp(all_sensors, stype));
        col = cfg(p).id + 1;
        avail_matrix(row, col)   = 1;
        primary_matrix(row, col) = cfg(p).sensors(s).primary;
    end
end

figure('Name','R&D: Sensor & Actuator Assignment per Phase','NumberTitle','off', ...
       'Position',[200 400 900 600]);

subplot(1,2,1);
imagesc(avail_matrix + primary_matrix);
colormap(subplot(1,2,1), [1 1 1; 0.7 0.85 1; 0.2 0.5 0.9]);
set(gca,'XTick',1:7,'XTickLabel',{'Ph0','Ph1','Ph2','Ph3','Ph4','Ph5','Ph6'});
set(gca,'YTick',1:numel(all_sensors),'YTickLabel', strrep(all_sensors,'_',' '));
xlabel('Phase'); title('Sensor Matrix (dark=primary, light=secondary)');
colorbar('Ticks',[0 1 2],'TickLabels',{'N/A','Secondary','Primary'});

% Actuator matrix
all_acts = {};
for p = 1:numel(cfg)
    for a = 1:numel(cfg(p).actuators)
        atype = cfg(p).actuators(a).type;
        if ~any(strcmp(all_acts, atype))
            all_acts{end+1} = atype; %#ok<AGROW>
        end
    end
end

act_matrix = zeros(numel(all_acts), 7);
for p = 1:numel(cfg)
    for a = 1:numel(cfg(p).actuators)
        atype = cfg(p).actuators(a).type;
        row = find(strcmp(all_acts, atype));
        col = cfg(p).id + 1;
        act_matrix(row, col) = 1;
    end
end

subplot(1,2,2);
imagesc(act_matrix);
colormap(subplot(1,2,2), [1 1 1; 0.2 0.7 0.3]);
set(gca,'XTick',1:7,'XTickLabel',{'Ph0','Ph1','Ph2','Ph3','Ph4','Ph5','Ph6'});
set(gca,'YTick',1:numel(all_acts),'YTickLabel', strrep(all_acts,'_',' '));
xlabel('Phase'); title('Actuator Matrix (green=active)');

sgtitle('Sensor & Actuator Assignment – Fehse (2003) R&D Phases', ...
        'FontWeight','bold','FontSize',12);

fprintf('[RD_PlotResults] All figures generated.\n');

end
