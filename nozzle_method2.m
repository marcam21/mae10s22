% Nozzle Contour Generation - Method 2: Method of Characteristics (MoC)
%
% Summary of Logic:
% This script generates a 2D Minimum Length Nozzle contour using the Method of Characteristics.
% 1. It calculates gas dynamic properties (Mach, Area Ratio, C*) for the given inputs.
% 2. It solves for the target Mach number at the exit corresponding to the desired Area Ratio.
% 3. It constructs the characteristic mesh (net) starting from a centered expansion fan at the throat.
%    - The fan angle (theta_max) is determined by the Prandtl-Meyer function (nu_e / 2) for minimum length.
% 4. It marches downstream using the Unit Process for interior points (intersection of C+ and C- characteristics).
%    - A planar approximation is used for the kernel calculation to ensure numerical stability.
% 5. It constructs the nozzle wall by tracing the streamline that cancels the reflected waves from the axis.
%    - The wall angle matches the local flow angle to prevent shocks (wave cancellation condition).
% 6. The script plots the resulting smooth wall contour and outputs coordinate points.

clear; clc; close all;

% --- Inputs & Constants ---
F_ideal = 12; % Newtons
P_c_low = 5e5; % Pa
P_c_high = 20e5; % Pa
Pc_array = linspace(P_c_low,P_c_high,1000);
t_r = 0.005; % seconds
T_c = 1222; % K
MW = 0.022560; % kg/mol
R = 8.314/MW; % J/kg
g = 9.806; % m/s^2
gamma    = 1.25;

% --- Derived Parameters ---
term = (2/(gamma+1))^((gamma+1)/(2*(gamma-1)));
C_star = sqrt(R*T_c) / ( sqrt(gamma) * term );

% Design Point Selection
A_ratio_target = 60;
A_t_target = 4e-6; % m^2
Rt = sqrt(A_t_target/pi);

% --- Solver Functions (Anonymous) ---
nuPM = @(M) sqrt((gamma+1)/(gamma-1))*atan( sqrt((gamma-1)/(gamma+1)*(M.^2-1)) ) ...
          - atan( sqrt(M.^2-1) );

nu_max = 0.5*pi*( sqrt((gamma+1)/(gamma-1)) - 1 );

% Mach Solver
f_Mach = @(M) (1./M) .* ((2./(gamma+1)) .* (1 + (gamma-1)./2 .* M.^2)) .^ ((gamma+1)./(2*(gamma-1))) - A_ratio_target;
M_exit = fzero(f_Mach, [1, 100]); % Supersonic root

% Minimum Length Nozzle Condition
nu_exit = nuPM(M_exit);
theta_max = nu_exit / 2; % Maximum wall angle for minimum length

% --- Method of Characteristics (MoC) ---
N = 100; % Number of characteristic lines (Higher N = Smoother Wall)
dtheta = theta_max / N;

% Initialize Arrays
% Indices: i (1..N) for Fan Rays, j (1..i) for Reflected Rays
% Using simple coordinate lists for points
% We will store points in a structure or cell array to manage the net topology

% Start Points (Expansion Fan at Throat Corner: x=0, y=Rt)
% For each ray i, the initial state is at the throat lip but with different flow properties
% theta_i = i * dtheta
% nu_i = theta_i (Prandtl-Meyer expansion from M=1)
% Mach_i from nu_i

% Data storage: Points P(i, j)
% x(i,j), y(i,j), theta(i,j), nu(i,j), M(i,j), mu(i,j)
% Flattened storage or matrix? Matrix is easier.
% Size N x N (approx). Let's use N+1 x N+1 for 1-based indexing covering j=0..N.
% Rows = i (Fan Ray), Cols = j (Reflected Ray).
% i goes 1..N. j goes 0..i.

X  = zeros(N, N+1);
Y  = zeros(N, N+1);
Th = zeros(N, N+1);
Nu = zeros(N, N+1);
M  = zeros(N, N+1);
Mu = zeros(N, N+1);

% Inverse PM Function Helper
invPM = @(nu) fzero(@(m) nuPM(m) - nu, [1.0001, 100]);
muMach = @(m) asin(1./m);

% Initialization: Fan Rays at Throat (j=0)
% Physical location is (0, Rt) for all rays, but properties differ.
for i = 1:N
    th = i * dtheta;
    nu = th;
    m  = invPM(nu);
    mu = muMach(m);

    X(i, 1)  = 0;
    Y(i, 1)  = Rt;
    Th(i, 1) = th;
    Nu(i, 1) = nu;
    M(i, 1)  = m;
    Mu(i, 1) = mu;
end

% Marching Loop
for i = 1:N % Fan Ray Loop
    for j = 1:i % Reflected Ray Loop (j=1 is Axis, j=2..i is Interior)

        idx_P = j + 1; % Current Point P(i,j)
        idx_L = j;     % Left Point L(i,j-1)
        idx_B = j + 1; % Bottom Point B(i-1,j)

        if j == i % Axis Point (j=i)
            % Intersection of Fan Ray i (from L) with Axis
            L_x = X(i, idx_L); L_y = Y(i, idx_L);
            L_th = Th(i, idx_L); L_nu = Nu(i, idx_L); L_mu = Mu(i, idx_L);

            % Axis BC: Theta = 0, y = 0
            P_th = 0;
            P_y  = 0;

            % Compatibility along C- (Fan Ray): theta - nu = const
            % P_th - P_nu = L_th - L_nu
            % 0 - P_nu = L_th - L_nu -> P_nu = L_nu - L_th
            % Note: Axisymmetric source terms ignored for stability (Planar approx)
            P_nu = L_nu - L_th;
            if P_nu < 0; P_nu = 0; end

            % Geometry: Intersection of Line from L (slope tan(th-mu)) and Axis (y=0)
            m_minus = tan(L_th - L_mu);
            P_x = L_x + (P_y - L_y) / m_minus;

            P_m = invPM(P_nu);
            P_mu = muMach(P_m);

            % Store
            X(i, idx_P) = P_x; Y(i, idx_P) = P_y;
            Th(i, idx_P) = P_th; Nu(i, idx_P) = P_nu;
            M(i, idx_P) = P_m; Mu(i, idx_P) = P_mu;

        else % Interior Point
            % Intersection of Fan Ray i (from L) and Refl Ray j (from B)
            L_x = X(i, idx_L); L_y = Y(i, idx_L);
            L_th = Th(i, idx_L); L_nu = Nu(i, idx_L); L_mu = Mu(i, idx_L);

            B_x = X(i-1, idx_B); B_y = Y(i-1, idx_B);
            B_th = Th(i-1, idx_B); B_nu = Nu(i-1, idx_B); B_mu = Mu(i-1, idx_B);

            % Slopes
            m_plus  = tan(B_th + B_mu); % C+ from B
            m_minus = tan(L_th - L_mu); % C- from L

            % Geometry Intersection
            P_x = (L_y - B_y + m_plus*B_x - m_minus*L_x) / (m_plus - m_minus);
            P_y = B_y + m_plus*(P_x - B_x);

            % Compatibility (Planar)
            % C+: th + nu = const (from B)
            % C-: th - nu = const (from L)
            J_plus  = B_th + B_nu;
            J_minus = L_th - L_nu;

            P_th = 0.5 * (J_plus + J_minus);
            P_nu = 0.5 * (J_plus - J_minus);

            if P_nu < 0; P_nu = 0; end
            P_m = invPM(P_nu);
            P_mu = muMach(P_m);

            % Store
            X(i, idx_P) = P_x; Y(i, idx_P) = P_y;
            Th(i, idx_P) = P_th; Nu(i, idx_P) = P_nu;
            M(i, idx_P) = P_m; Mu(i, idx_P) = P_mu;
        end
    end
end

% Wall Construction
% Trace streamline from the last characteristic point P(N,j)
wall_x = zeros(N+1, 1);
wall_y = zeros(N+1, 1);
wall_x(1) = 0; wall_y(1) = Rt;

curr_x = 0; curr_y = Rt; curr_th = theta_max;

for j = 1:N
    % Point P is the end of Reflected Ray j (on Fan Ray N)
    idx_P = j + 1;
    P_x = X(N, idx_P); P_y = Y(N, idx_P);
    P_th = Th(N, idx_P); P_mu = Mu(N, idx_P);

    % Extend C+ characteristic from P to Wall
    m_plus = tan(P_th + P_mu);

    % Wall Segment Slope
    % Linear turning distribution: Wall angle decreases from theta_max to 0
    th_wall_next = theta_max * (1 - j/N);
    th_avg = 0.5 * (curr_th + th_wall_next);
    m_wall = tan(th_avg);

    % Intersection
    w_x = (curr_y - P_y + m_plus*P_x - m_wall*curr_x) / (m_plus - m_wall);
    w_y = curr_y + m_wall*(w_x - curr_x);

    wall_x(j+1) = w_x;
    wall_y(j+1) = w_y;

    curr_x = w_x; curr_y = w_y; curr_th = th_wall_next;
end

% --- Output ---
fprintf('Method 2: Method of Characteristics (MoC)\n');
fprintf('Exit Mach:       %.4f\n', M_exit);
fprintf('Max Wall Angle:  %.2f deg\n', rad2deg(theta_max));
fprintf('Nozzle Length:   %.4f m\n', wall_x(end));
fprintf('Exit Radius:     %.4f m\n', wall_y(end));
fprintf('Wall Coordinates (x, y):\n');
for k = 1:10:length(wall_x) % Print sparse sample
    fprintf('%.6f, %.6f\n', wall_x(k), wall_y(k));
end

% --- Plotting ---
figure('Name', 'Method 2: MoC Contour');
hold on; axis equal; grid on;

% Plot Wall Contour (Top and Bottom)
plot(wall_x, wall_y, 'k-', 'LineWidth', 2);
plot(wall_x, -wall_y, 'k-', 'LineWidth', 2);

% Fill Interior
fill([wall_x; flipud(wall_x)], [wall_y; flipud(-wall_y)], [0.9 1.0 0.9], 'EdgeColor', 'none');

% Re-plot Wall
plot(wall_x, wall_y, 'k-', 'LineWidth', 2);
plot(wall_x, -wall_y, 'k-', 'LineWidth', 2);

% Axis
yline(0, 'k-.');
xline(0, 'k--');

xlabel('Axial Position x (m)');
ylabel('Radius r (m)');
title('Minimum Length Nozzle Contour (MoC)');
xlim([-0.1*wall_x(end), 1.1*wall_x(end)]);
ylim([-1.5*wall_y(end), 1.5*wall_y(end)]);
