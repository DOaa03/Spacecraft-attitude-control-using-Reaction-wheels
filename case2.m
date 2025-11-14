clear; clc; 

%% ===================== ORBIT PROPAGATION ===================== %%
mu = 398600; % Earth's gravitational parameter [km^3/s^2]
R0 = [3829.45; -888.41; 5459.13];    % Initial position [km]
V0 = [2.5396; 7.2434; -0.6063];      % Initial velocity [km/s]
y0 = [R0; V0];
tspan = [0 200];                    % Simulation time [s]

% Use ode45 with tighter tolerances
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t, y] = ode45(@(t, y) twoBodyODE(t, y, mu), tspan, y0, options);

r = y(:,1:3);
v = y(:,4:6);
N = length(t);

%% ===================== SPACECRAFT & WHEELS ===================== %%
J = diag([4.0, 3.0, 3.0]);     % Body inertia [kg·m²]
invJ = inv(J);
I = [0.05; 0.05; 0.05];        % Wheel inertia [kg·m²]
bw = [1e-3; 1e-3; 1e-3];       % Damping coefficients
S = eye(3);                    % Spin axes
q = [1; 0; 0; 0];              % Initial quaternion (scalar first)
omega = [0; 0; 0];             % Initial body rate
Omega = [0; 0; 0];             % Initial wheel speed

Kp = 0.5;                      % Proportional gain (كان 5.0)
Kd = 5.0;                      % Derivative gain (كان 10.0)

% Fixed time step for attitude and wheel dynamics
dt_fixed = 0.01; % Small fixed time step [s]
max_torque = 0.1; % Torque limit [Nm]


% % target case2
q_ref = [0.5,0.5,0.5,0.5];
omega_ref = [0;0;0];


%% ===================== STORAGE & INITIALIZATION ===================== %%
t_attitude = 0:dt_fixed:tspan(2);
N_attitude = length(t_attitude);

q_final = zeros(N_attitude, 4);
omega_final = zeros(N_attitude, 3);
Omega_final = zeros(N_attitude, 3);
tau_final = zeros(N_attitude, 3);
H_w_final = zeros(N_attitude, 3);
q_ref_final = zeros(N_attitude, 4);
q_err_final = zeros(N_attitude, 4);

q_final(1,:) = q';
omega_final(1,:) = omega';
Omega_final(1,:) = Omega';
q_ref_final(1,:) = q_ref ;
q_err_final(1,:) = [0.5,-0.5,-0.5,-0.5] ;

%% ===================== MAIN SIMULATION LOOP ===================== %%
i = 1; 
for k = 1:N-1 
    
  

   
    % دمج خطوات زمنية صغيرة وثابتة للتحكم والديناميكيات
    num_substeps = ceil((t(k+1) - t(k)) / dt_fixed);
    
    for sub = 1:num_substeps
        i = i + 1; 
        if i > N_attitude
            break;
        end
        
        % 1. Quaternion error
        q_err = quatmultiply(quatconj(q_ref), q');
        if q_err(1) < 0
            q_err = -q_err;
        end
        q_err = q_err'; 
        
        % 2. Angular velocity error
        omega_error = omega - omega_ref;
        % 3. PD Control torque (التصحيح: بدون حد إلغاء ذاتي + PD على الخطأ)
        % tau_controller = cross(omega, J*omega) - Kp * q_err(2:4) - Kd * omega_error; (النسخة القديمة)
        tau_controller = skew(omega)*J*omega  - Kp * q_err(2:4) - Kd * omega_error;
        
        % 4.
        tau_controller = min(max(tau_controller, -max_torque), max_torque);
        
        % 5.
        tau_RW = pinv(S) * (-tau_controller)
        
        tau_body = tau_controller;
        
        % 6. Wheel dynamics
        Omega_dot = (tau_RW - bw .* Omega) ./ I;
        Omega = Omega + Omega_dot * dt_fixed;
        
        % 7. Wheel angular momentum
        H_w = S * (I .* Omega);
        
        % 8. Body dynamics 
        IwOmegadot = S * (I .* Omega_dot);
        %J*omega_dot = tau_body - cross(omega, J*omega + H_w) - IwOmegadot
        omega_dot = invJ * (tau_body - cross(omega, J*omega + H_w) - IwOmegadot);
        
        omega = omega + omega_dot * dt_fixed;
        
        % 9. Quaternion kinematics
        capOmega = [ 0,         -omega(1), -omega(2), -omega(3);
            omega(1),  0,         omega(3), -omega(2);
            omega(2), -omega(3),  0,         omega(1);
            omega(3),  omega(2), -omega(1),  0];
        q_dot = 0.5 * capOmega * q;
        q = q + q_dot * dt_fixed
        
        % 10. Robust quaternion normalization
        q = q / norm(q);
        
        % Store results
        q_final(i,:) = q';
        q_ref_final(i,:) = q_ref';
        q_err_final(i,:) = q_err';
        omega_final(i,:) = omega';
        Omega_final(i,:) = Omega';
        tau_final(i,:) = tau_controller';
        H_w_final(i,:) = H_w';
    end
end

q_final = q_final(1:i-1,:);
q_ref_final = q_ref_final(1:i-1,:);
q_err_final = q_err_final(1:i-1,:);
omega_final = omega_final(1:i-1,:);
Omega_final = Omega_final(1:i-1,:);
tau_final = tau_final(1:i-1,:);
H_w_final = H_w_final(1:i-1,:);
t_attitude = t_attitude(1:i-1);

%% ===================== PLOTS ===================== %%

% current orientation and the refrence orientation
figure('Color','w');
subplot(3,1,1)
hold on;
plot(t_attitude,  q_final(:,1), 'k--','LineWidth', 1.8);
plot(t_attitude,  q_final(:,2), 'k-','LineWidth', 1.8);
plot(t_attitude,  q_final(:,3), 'k-.','LineWidth', 1.8);
plot(t_attitude,  q_final(:,4), 'k:','LineWidth', 1.8);
xlabel('Time [s]'); ylabel('q');
legend('q_0','q_1','q_2','q_3','q_{0,des}','q_{1,des}','q_{2,des}','q_{3,des}', 'Location', 'best');
title('Quaternion Evolution');
grid on;

subplot(3,1,2)

hold on
plot(t_attitude,  q_err_final(:,1), 'k--','LineWidth', 1.8);
plot(t_attitude,  q_err_final(:,2), 'k-','LineWidth', 1.8);
plot(t_attitude,  q_err_final(:,3), 'k-.','LineWidth', 1.8);
plot(t_attitude,  q_err_final(:,4), 'k:','LineWidth', 1.8);
xlabel('Time [s]');
legend('q_{err0}','q_{err1}','q_{err2}','q_{err3}', 'Location', 'best');
title('Quaternion Error ');
grid on;


figure('Color','w');
subplot(3,1,1);
plot(t_attitude, omega_final(:,1),'k','LineWidth',1.8);
xlabel('Time [s]'); ylabel('\omega_x [rad/s]');
title('Body Angular Velocity - \omega_x');
grid on;

subplot(3,1,2);
plot(t_attitude, omega_final(:,2),'k','LineWidth',1.8);
xlabel('Time [s]'); ylabel('\omega_y [rad/s]');
title('Body Angular Velocity - \omega_y');
grid on;

subplot(3,1,3);
plot(t_attitude, omega_final(:,3),'k','LineWidth',1.8);
xlabel('Time [s]'); ylabel('\omega_z [rad/s]');
title('Body Angular Velocity - \omega_z');
grid on;

sgtitle('Spacecraft Body Angular Velocities (No Control)');

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

sgtitle('Reaction Wheel Angular Velocities (No Control)');

% -----------control torque -----------
figure('Color','w');
subplot(3,1,1)
plot(t_attitude, vecnorm(tau_final,2,2),'k','LineWidth',1.8);
xlabel('Time [s]'); ylabel('|\tau_{controller}| [N·m]');
title('Required control torque');
grid on;

%------------reaction wheel angular momentom-----------------
subplot(3,1,2)
plot(t_attitude, vecnorm(H_w_final,2,2),'k','LineWidth',1.8);
xlabel('Time [s]'); ylabel('|H_{w}| [N·m.s]');
title("Wheels angular momentum");
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

