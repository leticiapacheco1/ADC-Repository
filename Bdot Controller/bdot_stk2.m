%% DARLA 02 Bdot Sim
%  Runs Euler+quaternion dynamics in MATLAB (RK4)
%  Uses igrfmagm() if available, falls back to dipole model
%  Detects detumble (<5 deg/s on every axis) and confirms when held for 1 full orbit

clc; clear; close all;

%% mission parameters 

STK_HOST  = 'localhost';
STK_PORT  = 5001;
scenName  = 'DARLA02_current';
satName   = 'DARLA02';

startEpoch = datetime(2026, 6, 5, 0, 0, 0);
dyear = year(startEpoch) + day(startEpoch,'dayofyear')/365;

simDays  = 4;
dt       = 0.5;
TOTAL_N  = round(simDays*86400/dt);

Ixx = 0.024;  Iyy = 0.023;  Izz = 0.006;
I     = diag([Ixx, Iyy, Izz]);
I_inv = inv(I);

% mtq sizing — 26 AWG, 6 layers, 3U face
% typical 3U face area ≈ 0.01 × 0.1 m² = 1e-3 m² per layer pair
% m = N·I·A  -> (0.3 A) × (500 turns) × (1e-3 m² × 6-layer)
m_max = 0.15;                       % [A·m²] per axis


tipoff_degs = 11.0;
omega0 = deg2rad([tipoff_degs; -tipoff_degs*0.9; tipoff_degs*0.95]);

% ADC requirement threshold
thresh_req  = deg2rad(5.0);         % [rad/s]
thresh_goal = deg2rad(0.75);        % [rad/s]

% orbit
mu     = 398600.4418;               % [km³/s²]
a_orb  = 6878.137;                  % [km] 500 km alt approximatelly 
n_orb  = sqrt(mu / a_orb^3);       % [rad/s] mean motion
T_orb  = 2*pi / n_orb;             % [s] orbital period
inc    = deg2rad(63.4);             % inclination

% b-dot gain
% k = 2·n·(1+sin(i))·λ_max(I) / B_avg²  × tuning
B_avg  = 3.0e-5;                    % [T] average LEO field
k_bdot = 2 * n_orb * (1+sin(inc)) * max([Ixx Iyy Izz]) / B_avg^2;
k_bdot = k_bdot * 0.30; 

fprintf('DARLA02 B-DOT DETUMBLING\n');
fprintf('Tip-off:     %.1f deg/s per axis\n', tipoff_degs);
fprintf('m_max:       %.3f A·m² per axis\n', m_max);
fprintf('k_bdot:      %.4e\n', k_bdot);
fprintf('Req:         < %.1f deg/s  |  Goal: < %.2f deg/s\n', ...
    rad2deg(thresh_req), rad2deg(thresh_goal));
fprintf('Orb period:  %.2f min\n', T_orb/60);

%% stk
fprintf('\nConnecting to STK12 (%s:%d)...\n', STK_HOST, STK_PORT);
tcp = tcpclient(STK_HOST, STK_PORT, 'Timeout', 10);
pause(0.5);
if tcp.NumBytesAvailable > 0
    read(tcp, tcp.NumBytesAvailable);   % flush banner
end
fprintf('Connected.\n\n');

%% log
time_log  = zeros(TOTAL_N, 1);
omega_log = zeros(TOTAL_N, 3);      % [deg/s]

%% animated plot
hFig = figure('Name','DARLA02 B-dot Detumbling — Live','Color','w', ...
              'Position',[60 60 920 520]);

subplot(2,1,1);   % main angular velocity plot
hWx = animatedline('Color','b','LineWidth',1.2,'DisplayName','\omega_x');
hWy = animatedline('Color','r','LineWidth',1.2,'DisplayName','\omega_y');
hWz = animatedline('Color',[0.85 0.65 0],'LineWidth',1.2,'DisplayName','\omega_z');
hReq = yline(rad2deg(thresh_req),'--','Color',[0.2 0.6 0.2],'LineWidth',1.2, ...
    'Label','+5 deg/s req.');
yline(-rad2deg(thresh_req),'--','Color',[0.2 0.6 0.2],'LineWidth',1.2, ...
    'Label','-5 deg/s req.');
xlabel('Time [min]','FontSize',10);
ylabel('Angular rate [deg/s]','FontSize',10);
title('DARLA02 — B-dot Detumbling (live)','FontSize',12,'FontWeight','bold');
legend([hWx hWy hWz],'Location','northeast','FontSize',10);
ylim([-15 15]);  grid on;  box on;

subplot(2,1,2);   % norm
hNorm = animatedline('Color','k','LineWidth',1.2);
yline(rad2deg(thresh_req),'--','Color',[0.2 0.6 0.2],'LineWidth',1.2);
xlabel('Time [min]','FontSize',10);
ylabel('|\omega| [deg/s]','FontSize',10);
title('Angular rate magnitude','FontSize',11);
ylim([0 25]); grid on; box on;

%% initial State 
omega  = omega0;
q      = [0; 0; 0; 1];
B_prev = zeros(3,1);

% low-pass filter for Bdot
tau_lp   = 20;                      % [s]
alpha_lp = dt / (tau_lp + dt);
Bdot_filt = zeros(3,1);

% detumble state
detumbled    = false;
confirmed    = false;
detumble_t   = NaN;
confirm_t    = NaN;
req_hold_cnt = 0;                   % consecutive steps below threshold
req_hold_need = round(T_orb / dt); % must hold for 1 orbit

%% main loop 
fprintf('Running simulation... (STK animates live)\n');
plot_every = 20;     % update plot every N steps
stk_every  = 2;      % send to STK every N steps (don't flood TCP)

for k = 1:TOTAL_N

    t_now = (k-1) * dt;            % [s]
    t_min = t_now / 60;            % [min] for logging

    % orbital position (lat/lon for IGRF)
    u   = n_orb * t_now;           % argument of latitude [rad]
    lat = rad2deg(asin(sin(inc) * sin(u)));
    lon = rad2deg(atan2(cos(inc)*sin(u), cos(u)));
    alt_km = 500;

    % B-field: igrfmagm, dipole fallback
    try
        [Bn, Be, Bd] = igrfmagm(lat, lon, alt_km, dyear);
        if any(isnan([Bn Be Bd])), error('nan'); end
    catch
        [Bn, Be, Bd] = dipole_B(lat, alt_km);
    end

    slat=sind(lat); clat=cosd(lat); slon=sind(lon); clon=cosd(lon);
    R_ned2ecef = [-slat*clon, -slon, -clat*clon;
                  -slat*slon,  clon, -clat*slon;
                   clat,       0,    -slat];
    B_ecef = R_ned2ecef * [Bn; Be; Bd] * 1e-9;  % [T]

    % quaternion
    qw=q(4); qx=q(1); qy=q(2); qz=q(3);
    R_b2i = [1-2*(qy^2+qz^2),   2*(qx*qy-qz*qw), 2*(qx*qz+qy*qw);
             2*(qx*qy+qz*qw),   1-2*(qx^2+qz^2), 2*(qy*qz-qx*qw);
             2*(qx*qz-qy*qw),   2*(qy*qz+qx*qw), 1-2*(qx^2+qy^2)];
    B_body = R_b2i' * B_ecef;

    if any(isnan(B_body)), B_body = B_prev; end

    % Bdot w/ low pass filter (to prevent spikes)
    Bdot_raw  = (B_body - B_prev) / dt;
    Bdot_filt = alpha_lp * Bdot_raw + (1 - alpha_lp) * Bdot_filt;

    m = -k_bdot * Bdot_filt;   % Bdot control law 
    m = max(-m_max, min(m_max, m));

    tau = cross(m, B_body); % magnetic torque
    omega = rk4_omega(omega, tau, I, I_inv, dt); % RK4 angular velocity

    q = quat_rk4(q, omega, dt);

    
    od = rad2deg(omega);
    time_log(k)    = t_min;
    omega_log(k,:) = od';
    B_prev = B_body;

    all_below_req = all(abs(omega) < thresh_req);
    if ~detumbled
        if all_below_req
            detumbled  = true;
            detumble_t = t_now;
            req_hold_cnt = 1;
            fprintf('  >> Below 5 deg/s at %.2f orbits = %.1f min\n', ...
                t_now/T_orb, t_min);
        end
    else
        if all_below_req
            req_hold_cnt = req_hold_cnt + 1;
            if ~confirmed && req_hold_cnt >= req_hold_need
                confirmed = true;
                confirm_t = t_now;
                fprintf('  >> DETUMBLE CONFIRMED at %.2f orbits = %.1f hrs\n', ...
                    t_now/T_orb, t_now/3600);
            end
        else
            % re-tumbled
            detumbled    = false;
            req_hold_cnt = 0;
            detumble_t   = NaN;
        end
    end

    % send attitude to stk
    if mod(k, stk_every) == 0
        curTime = startEpoch + seconds(t_now);
        timeStr = stk_time_str(curTime);

        % per-step quaternion attitude command
        cmd = sprintf('SetAttitude */Satellite/%s Quaternion "%s" %.8f %.8f %.8f %.8f', ...
            satName, timeStr, qw, qx, qy, qz);
        write(tcp, [uint8(cmd), 13, 10]);

        cmd2 = sprintf('Animate */Scenario/%s SetTime "%s"', scenName, timeStr);
        write(tcp, [uint8(cmd2), 13, 10]);
    end

    
    if mod(k, plot_every) == 0
        addpoints(hWx, t_min, od(1));
        addpoints(hWy, t_min, od(2));
        addpoints(hWz, t_min, od(3));
        addpoints(hNorm, t_min, norm(od));
        drawnow limitrate;
    end

    if confirmed && all(abs(omega) < thresh_goal)
        fprintf('  >> Goal (%.2f deg/s) reached at %.2f orbits\n', ...
            rad2deg(thresh_goal), t_now/T_orb);
        TOTAL_N = k;   % trim log
        break;
    end
end

%% final plot
t_log  = time_log(1:TOTAL_N);
wx_log = omega_log(1:TOTAL_N,1);
wy_log = omega_log(1:TOTAL_N,2);
wz_log = omega_log(1:TOTAL_N,3);

t_orb_log = t_log * 60 / T_orb;
ol = sprintf('(1 orbit = %.2f hrs; 1 day = %.2f orbits)', T_orb/3600, 86400/T_orb);

figure('Name','DARLA02 — Final Detumbling Result','Color','w','Position',[100 100 1000 540]);
plot(t_orb_log, wx_log,'b','LineWidth',1.2); hold on;
plot(t_orb_log, wy_log,'r','LineWidth',1.2);
plot(t_orb_log, wz_log,'Color',[0.85 0.65 0],'LineWidth',1.2);
yline( rad2deg(thresh_req),'--','Color',[0.2 0.6 0.2],'LineWidth',1.4, ...
    'Label','+5 deg/s req.','FontSize',9);
yline(-rad2deg(thresh_req),'--','Color',[0.2 0.6 0.2],'LineWidth',1.4, ...
    'Label','-5 deg/s req.','FontSize',9);
if ~isnan(detumble_t)
    xline(detumble_t/T_orb,'--k','LineWidth',1.2,'Alpha',0.5,'Label','Detumble');
end
if ~isnan(confirm_t)
    xline(confirm_t/T_orb,'--','Color',[0.2 0.6 0.2],'LineWidth',1.2, ...
        'Alpha',0.6,'Label','Confirmed');
    text(max(t_orb_log)*0.32, -11, ...
        sprintf(['\\bfDetumbled: %.1f orbits = %.1f hrs\n' ...
                 '\\bfConfirmed: %.1f orbits\n' ...
                 '\\bfThreshold: ±5 deg/s'], ...
        detumble_t/T_orb, detumble_t/3600, confirm_t/T_orb), ...
        'FontSize',9.5,'BackgroundColor','w','EdgeColor','k');
end
xlabel(['Number of orbits  ' ol],'FontSize',10);
ylabel('Angular rates [deg/s]','FontSize',11);
title(sprintf('DARLA02 B-dot Detumbling  |  Tip-off %.0f deg/s  |  m_{max}=%.2f A·m²', ...
    tipoff_degs, m_max),'FontSize',12,'FontWeight','bold');
legend('\omega_x','\omega_y','\omega_z','Location','northeast','FontSize',11);
ylim([-15 15]); grid on; box on;

saveas(gcf, 'DARLA02_detumbling_result.png');

%% excel
xls_every = max(1, round(30/dt));
xi = 1:xls_every:TOTAL_N;

Tout = table( ...
    time_log(xi),   ...
    omega_log(xi,1),...
    omega_log(xi,2),...
    omega_log(xi,3),...
    'VariableNames', {'Time_min','omega_x_degs','omega_y_degs','omega_z_degs'});
writetable(Tout, 'DARLA02_angular_velocity.xlsx', 'Sheet','Angular_Velocity');

detStr = 'Not reached';
confStr = 'Not reached';
if ~isnan(detumble_t)
    detStr = sprintf('%.2f orbits = %.1f hrs', detumble_t/T_orb, detumble_t/3600);
end
if ~isnan(confirm_t)
    confStr = sprintf('%.2f orbits', confirm_t/T_orb);
end

S = {'Parameter','Value'; ...
     'Satellite', satName; ...
     'Epoch', stk_epoch_str(startEpoch); ...
     'Tip-off [deg/s]', tipoff_degs; ...
     'Ixx [kg·m²]', Ixx; ...
     'Iyy [kg·m²]', Iyy; ...
     'Izz [kg·m²]', Izz; ...
     'm_max per axis [A·m²]', m_max; ...
     'B-dot gain k', k_bdot; ...
     'Threshold [deg/s]', 5.0; ...
     'Orbital period [min]', T_orb/60; ...
     'Detumbled at', detStr; ...
     'Confirmed stable at', confStr};
writetable(cell2table(S(2:end,:),'VariableNames',S(1,:)), ...
    'DARLA02_angular_velocity.xlsx','Sheet','Summary');

fprintf('\nExcel saved: DARLA02_angular_velocity.xlsx\n');

%% summary
fprintf('\n');
fprintf('DARLA02 DETUMBLING REPORT \n');
fprintf('Inertia:       Ixx=%.3f  Iyy=%.3f  Izz=%.3f kg·m²\n', Ixx,Iyy,Izz);
fprintf('MTQ m_max:     %.3f A·m² per axis (26 AWG, 6 layers)\n', m_max);
fprintf('Initial tumble: %.1f deg/s (tip-off)\n', tipoff_degs);
fprintf('ADC requirement: < 5 deg/s on all axes\n');

if ~isnan(detumble_t)
    fprintf('Detumbled at:  %.1f hrs  (%.1f orbits = %.2f days)\n', ...
        detumble_t/3600, detumble_t/T_orb, detumble_t/86400);
    fprintf('Confirmed at:  %.1f orbits  (stable for 1 full orbit)\n', ...
        confirm_t/T_orb);
    fprintf('RESULT: Depending on a %.0f deg/s tip-off rate,\n', tipoff_degs);
    fprintf('         DARLA02 will detumble in %.1f hours\n', detumble_t/3600);
    fprintf('         (%.1f days / %.0f orbits).\n', detumble_t/86400, detumble_t/T_orb);
    fprintf('         This meets the <5 deg/s ADC requirement.\n');
else
    fprintf(' RESULT: Detumble NOT achieved in %.0f days.\n', simDays);
    fprintf('         Increase m_max or check MTQ sizing.\n');
end

%% local functions

function w = rk4_omega(omega, tau, I, I_inv, dt)
    f  = @(w) I_inv * (tau - cross(w, I*w));
    k1 = f(omega);
    k2 = f(omega + 0.5*dt*k1);
    k3 = f(omega + 0.5*dt*k2);
    k4 = f(omega + dt*k3);
    w  = omega + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end

function q_new = quat_rk4(q, omega, dt)
    wx=omega(1); wy=omega(2); wz=omega(3);
    Om = [ 0,  wz, -wy, wx;
          -wz,  0,  wx, wy;
           wy, -wx,  0, wz;
          -wx, -wy, -wz, 0];
    f  = @(qq) 0.5*Om*qq;
    k1 = f(q);
    k2 = f(q + 0.5*dt*k1);
    k3 = f(q + 0.5*dt*k2);
    k4 = f(q + dt*k3);
    q_new = q + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    q_new = q_new / norm(q_new);
end

function [Bn, Be, Bd] = dipole_B(lat, alt_km)
    Re   = 6371; r = Re + alt_km; B0 = 3.12e4;
    latr = deg2rad(lat);
    Br   = -2*B0*(Re/r)^3*sin(latr);
    Bt   =    B0*(Re/r)^3*cos(latr);
    Bn   = -Bt;  Be = 0;  Bd = -Br;
end

function s = stk_time_str(dt_obj)
    mn = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
    s  = sprintf('%02d %s %04d %02d:%02d:%06.3f', ...
        day(dt_obj), mn{month(dt_obj)}, year(dt_obj), ...
        hour(dt_obj), minute(dt_obj), second(dt_obj));
end

function s = stk_epoch_str(dt_obj)
    s = stk_time_str(dt_obj);
end
