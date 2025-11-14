% هذا الكود يحاكي توجيه مركبة فضائية باستخدام ثلاث عجلات تفاعلية (Reaction Wheels)
% باستخدام قانون تحكم يعتمد على Lyapunov/PD لتتبع محور LVLH المداري
% النسخة المعدلة: تشمل النموذج النسبي للعجلات + حساب زخم العجلات والكلي
% بالإضافة إلى حساب وعرض السرعة الزاوية الحقيقية في إطار ECI

clear; clc;

%% ===================== ORBIT PROPAGATION ===================== %%
mu = 398600; % Earth's gravitational parameter [km^3/s^2]
R0 = [3829.45; -888.41; 5459.13];   % Initial position [km]
V0 = [2.5396; 7.2434; -0.6063];     % Initial velocity [km/s]
y0 = [R0; V0];
tspan = [0 200];                    % Simulation time [s]

% Use ode45 with tighter tolerances
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t, y] = ode45(@(t, y) twoBodyODE(t, y, mu), tspan, y0, options);

r = y(:,1:3);
v = y(:,4:6);
N = length(t);

%% ===================== SPACECRAFT & WHEELS ===================== %%
J = diag([4.0, 5.0, 3.0]);      % Body inertia [kg·m²]
invJ = inv(J);
I = [0.05; 0.05; 0.05];         % Wheel inertia [kg·m²]
bw = [1e-3; 1e-3; 1e-3];        % Wheel damping coefficients
S = eye(3);                     % Wheel spin axes (orthogonal)
q = [1; 0; 0; 0];               % Initial quaternion (scalar first)
omega = [0; 0; 0];              % Initial body angular velocity [rad/s]
Omega_rel = [0; 0; 0];          % Initial wheel relative speeds [rad/s]

% PD control gains (tuned for stable LVLH tracking)
Kp = 0.5;
Kd = 5.0;

% Time step for attitude simulation
dt_fixed = 0.01;
max_torque = 1.0; % Torque saturation [Nm]

%% ===================== STORAGE & INITIALIZATION ===================== %%
t_attitude = 0:dt_fixed:tspan(2);
N_attitude = length(t_attitude);

q_final = zeros(N_attitude, 4);
omega_final = zeros(N_attitude, 3);
omega_des_body_final = zeros(N_attitude, 3);
Omega_final = zeros(N_attitude, 3);
tau_final = zeros(N_attitude, 3);
q_ref_final = zeros(N_attitude, 4);
q_err_final = zeros(N_attitude, 4);
H_w_final = zeros(N_attitude, 3);
H_total_final = zeros(N_attitude, 3);

% Arrays for ECI comparisons
omega_ECI_final = zeros(N_attitude, 3);
omega_des_ECI_final = zeros(N_attitude, 3);

% initial store
q_final(1,:) = q';
omega_final(1,:) = omega';
Omega_final(1,:) = Omega_rel';

% compute and store initial LVLH reference and error (t = 0)
r_vec0 = r(1,:)';
v_vec0 = v(1,:)';
r_norm0 = norm(r_vec0);
z_LVLH0 = -r_vec0 / r_norm0;
h_vec0 = cross(r_vec0, v_vec0);
h_norm0 = norm(h_vec0);
y_LVLH0 = -h_vec0 / h_norm0;
x_LVLH0 = cross(y_LVLH0, z_LVLH0);
DCM_ECI2LVLH0 = [x_LVLH0, y_LVLH0, z_LVLH0];

q_ref0 = dcm2quat(DCM_ECI2LVLH0);
if dot(q_ref0, q_final(1,:)) < 0
    q_ref0 = -q_ref0;
end
q_ref0 = q_ref0';
q_ref_final(1,:) = q_ref0';

q_err0 = quatmultiply(quatconj(q_ref0'), q');
if q_err0(1) < 0
    q_err0 = -q_err0;
end
q_err0 = q_err0';
q_err_final(1,:) = q_err0';

omega_orbit0 = h_norm0 / r_norm0^2;
omega_des_LVLH0 = [0; -omega_orbit0; 0];
omega_des_body0 = DCM_ECI2LVLH0' * omega_des_LVLH0;
omega_des_body_final(1,:) = omega_des_body0';

% store initial ECI quantities (approx)
omega_ECI_final(1,:) = (DCM_ECI2LVLH0 * omega)';          % approx mapping of omega to ECI basis
omega_des_ECI_final(1,:) = (DCM_ECI2LVLH0 * omega_des_LVLH0)';

%% ===================== MAIN SIMULATION LOOP ===================== %%
i = 1;
for k = 1:N-1

    % تحديث إطار LVLH المرجعي بناءً على المدار الحالي
    r_vec = r(k,:)';
    v_vec = v(k,:)';
    r_norm = norm(r_vec);

    z_LVLH = -r_vec / r_norm;
    h_vec = cross(r_vec, v_vec);
    h_norm = norm(h_vec);
    y_LVLH = -h_vec / h_norm;
    x_LVLH = cross(y_LVLH, z_LVLH);
    DCM_ECI2LVLH = [x_LVLH, y_LVLH, z_LVLH];

    omega_orbit = h_norm / r_norm^2;
    omega_des_LVLH = [0; -omega_orbit; 0]; % reference rotation rate (LVLH)

    % Quaternion reference
    q_ref = dcm2quat(DCM_ECI2LVLH);
    if dot(q_ref, q_final(i,:)) < 0
        q_ref = -q_ref;
    end
    q_ref = q_ref';

    omega_des_body = DCM_ECI2LVLH' * omega_des_LVLH;

    % عدد الخطوات الصغيرة لكل فترة مدارية
    num_substeps = ceil((t(k+1) - t(k)) / dt_fixed);

    for sub = 1:num_substeps
        i = i + 1;
        if i > N_attitude
            break;
        end

        %% 1. Quaternion error
        q_err = quatmultiply(quatconj(q_ref'), q');
        if q_err(1) < 0
            q_err = -q_err;
        end
        q_err = q_err';

        %% 2. Angular velocity error
        omega_error = omega - omega_des_body;

        %% 3. PD Control torque
        tau_controller = skew(omega)*J*omega  - Kp * q_err(2:4) - Kd * omega_error;
        tau_controller = min(max(tau_controller, -max_torque), max_torque);

        %% 4. Torque on wheels (input torque)
        tau_RW = pinv(S) * (-tau_controller);

        %% 5. Compute wheel angular momentum (using relative model)
        Omega_abs = Omega_rel + omega;
        H_w = S * (I .* Omega_abs);

        %% 6. Body dynamics
        omega_dot = invJ * (tau_controller - cross(omega, J*omega + H_w));

        %% 7. Reaction wheel dynamics (relative model)
        Omega_rel_dot = (tau_RW - bw .* Omega_rel - (I .* omega_dot)) ./ I;
        Omega_rel = Omega_rel + Omega_rel_dot * dt_fixed;

        %% 8. Update absolute wheel speed and body rate
        Omega_abs = Omega_rel + omega;
        omega = omega + omega_dot * dt_fixed;

        %% 9. Quaternion kinematics
        capOmega = [ 0,         -omega(1), -omega(2), -omega(3);
                     omega(1),  0,         omega(3), -omega(2);
                     omega(2), -omega(3),  0,         omega(1);
                     omega(3),  omega(2), -omega(1),  0];
        q_dot = 0.5 * capOmega * q;
        q = q + q_dot * dt_fixed;
        q = q / norm(q);

        %% 10. Compute total angular momentum
        H_total = J*omega + H_w;

        %% 11. Store results
        q_final(i,:) = q';
        q_ref_final(i,:) = q_ref';
        q_err_final(i,:) = q_err';
        omega_final(i,:) = omega';
        omega_des_body_final(i,:) = omega_des_body';
        Omega_final(i,:) = Omega_abs';
        tau_final(i,:) = tau_controller';
        H_w_final(i,:) = H_w';
        H_total_final(i,:) = H_total';

        % ---- store ECI-referenced angular velocities for comparison ----
        % approximate mapping of body angular velocity into ECI basis using LVLH axes
        omega_ECI = DCM_ECI2LVLH * omega;               % approx mapping to ECI basis
        omega_des_ECI = DCM_ECI2LVLH * omega_des_LVLH;  % desired (LVLH -> ECI)

        omega_ECI_final(i,:) = omega_ECI';
        omega_des_ECI_final(i,:) = omega_des_ECI';
    end
end

% قص النتائج
q_final = q_final(1:i-1,:);
q_ref_final = q_ref_final(1:i-1,:);
q_err_final = q_err_final(1:i-1,:);
omega_final = omega_final(1:i-1,:);
omega_des_body_final = omega_des_body_final(1:i-1,:);
Omega_final = Omega_final(1:i-1,:);
tau_final = tau_final(1:i-1,:);
H_w_final = H_w_final(1:i-1,:);
H_total_final = H_total_final(1:i-1,:);
omega_ECI_final = omega_ECI_final(1:i-1,:);
omega_des_ECI_final = omega_des_ECI_final(1:i-1,:);
t_attitude = t_attitude(1:i-1);

%% ===================== PLOTS ===================== %%

% Quaternion LVLH ref and actual
figure('Color','w');
subplot(3,1,1)
hold on;
plot(t_attitude,  q_ref_final(:,1), 'k--','LineWidth', 1.2);
plot(t_attitude,  q_ref_final(:,2), 'k-','LineWidth', 1.2);
plot(t_attitude,  q_ref_final(:,3), 'k-.','LineWidth', 1.2);
plot(t_attitude,  q_ref_final(:,4), 'k:','LineWidth', 1.2);
xlabel('Time [s]'); ylabel('q (ref)');
title('Quaternion of LVLH through the orbit');
grid on;

subplot(3,1,2)
hold on;
plot(t_attitude,  q_final(:,1), 'k--','LineWidth', 1.8);
plot(t_attitude,  q_final(:,2), 'k-','LineWidth', 1.8);
plot(t_attitude,  q_final(:,3), 'k-.','LineWidth', 1.8);
plot(t_attitude,  q_final(:,4), 'k:','LineWidth', 1.8);
xlabel('Time [s]'); ylabel('q (body)');
title('Quaternion Evolution');
grid on;

subplot(3,1,3)
hold on
plot(t_attitude,  q_err_final(:,1), 'k','LineWidth', 1.8);
plot(t_attitude,  q_err_final(:,2), 'k','LineWidth', 1.8);
plot(t_attitude,  q_err_final(:,3), 'k','LineWidth', 1.8);
plot(t_attitude,  q_err_final(:,4), 'k','LineWidth', 1.8);
xlabel('Time [s]'); ylabel('q_{err}');
title('Quaternion Error');
grid on;

% Body angular velocities (body frame)
figure('Color','w');
subplot(3,1,1);
plot(t_attitude, omega_final(:,1),'k','LineWidth',1.8);
xlabel('Time [s]'); ylabel('\omega_x [rad/s]');
title('Body Angular Velocity - \omega_x (body frame)');
grid on;

subplot(3,1,2);
plot(t_attitude, omega_final(:,2),'k','LineWidth',1.8);
xlabel('Time [s]'); ylabel('\omega_y [rad/s]');
title('Body Angular Velocity - \omega_y (body frame)');
grid on;

subplot(3,1,3);
plot(t_attitude, omega_final(:,3),'k','LineWidth',1.8);
xlabel('Time [s]'); ylabel('\omega_z [rad/s]');
title('Body Angular Velocity - \omega_z (body frame)');
grid on;

sgtitle('Spacecraft Body Angular Velocities');

% Desired angular velocity in body frame
figure('Color','w');
subplot(3,1,1);
plot(t_attitude, omega_des_body_final(:,1),'r' ,'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('\omega_{des,x} [rad/s]');
title('Desired \omega_x in Body Frame');
grid on;

subplot(3,1,2);
plot(t_attitude, omega_des_body_final(:,2),'r' ,'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('\omega_{des,y} [rad/s]');
title('Desired \omega_y in Body Frame');
grid on;

subplot(3,1,3);
plot(t_attitude, omega_des_body_final(:,3),'r' ,'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('\omega_{des,z} [rad/s]');
title('Desired \omega_z in Body Frame');
grid on;

sgtitle('Desired Angular Velocity (Body Frame)');

% ---- ECI comparison: actual vs desired (ECI-referenced) ----
figure('Color','w');
subplot(3,1,1);
plot(t_attitude, omega_ECI_final(:,1),'k','LineWidth',1.2); hold on;
plot(t_attitude, omega_des_ECI_final(:,1),'r--','LineWidth',1.2);
xlabel('Time [s]'); ylabel('\omega_x [rad/s]'); legend('Actual (ECI approx)','Desired (ECI)');
title('Angular Velocity about X-axis (ECI basis)');
grid on;

subplot(3,1,2);
plot(t_attitude, omega_ECI_final(:,2),'k','LineWidth',1.2); hold on;
plot(t_attitude, omega_des_ECI_final(:,2),'r--','LineWidth',1.2);
xlabel('Time [s]'); ylabel('\omega_y [rad/s]'); legend('Actual (ECI approx)','Desired (ECI)');
title('Angular Velocity about Y-axis (ECI basis)');
grid on;

subplot(3,1,3);
plot(t_attitude, omega_ECI_final(:,3),'k','LineWidth',1.2); hold on;
plot(t_attitude, omega_des_ECI_final(:,3),'r--','LineWidth',1.2);
xlabel('Time [s]'); ylabel('\omega_z [rad/s]'); legend('Actual (ECI approx)','Desired (ECI)');
title('Angular Velocity about Z-axis (ECI basis)');
grid on;

sgtitle('Spacecraft True vs Desired Angular Velocity in ECI Basis (approx)');

% ----------- Reaction Wheel Velocities -----------
figure('Color','w');
subplot(3,1,1);
plot(t_attitude, Omega_final(:,1),'k','LineWidth',1.8);
xlabel('Time [s]'); ylabel('\Omega_x [rad/s]');
title('Reaction Wheel Speed - \Omega_x');
grid on;

subplot(3,1,2);
plot(t_attitude, Omega_final(:,2),'k','LineWidth',1.8);
xlabel('Time [s]'); ylabel('\Omega_y [rad/s]');
title('Reaction Wheel Speed - \Omega_y');
grid on;

subplot(3,1,3);
plot(t_attitude, Omega_final(:,3),'k','LineWidth',1.8);
xlabel('Time [s]'); ylabel('\Omega_z [rad/s]');
title('Reaction Wheel Speed - \Omega_z');
grid on;

sgtitle('Reaction Wheel Angular Velocities');

% -----------control torque -----------
figure('Color','w');
subplot(2,1,1)
plot(t_attitude, vecnorm(tau_final,2,2),'k','LineWidth',1.8);
xlabel('Time [s]'); ylabel('|\tau_{controller}| [N·m]');
title('Required control torque');
grid on;

%------------reaction wheel angular momentum-----------------
subplot(2,1,2)
plot(t_attitude, vecnorm(H_w_final,2,2),'k','LineWidth',1.8);
xlabel('Time [s]'); ylabel('|H_{w}| [N·m.s]');
title('Wheels angular momentum');
grid on;

%% ===================== ODE FUNCTION ===================== %%
function dydt = twoBodyODE(~, y, mu)
    x = y(1); y_pos = y(2); z = y(3);
    vx = y(4); vy = y(5); vz = y(6);
    r = sqrt(x^2 + y_pos^2 + z^2);
    if r < 1e-6
        r = 1e-6;
    end
    ax = -mu * x / r^3;
    ay = -mu * y_pos / r^3;
    az = -mu * z / r^3;
    dydt = [vx; vy; vz; ax; ay; az];
end

function S = skew(v)
S = [  0    -v(3)   v(2);
      v(3)    0    -v(1);
     -v(2)   v(1)    0  ];
end
