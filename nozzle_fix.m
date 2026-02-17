clear; clc; close all;

% Inputs
F_ideal = 12; % Newtons
P_c_low = 5e5; % Pa
P_c_high = 20e5; % Pa
Pc = linspace(P_c_low,P_c_high,1000);
t_r = 0.005; % seconds

% Constants
T_c = 1222; % K
MW = 0.022560; % kg/mol
R = 8.314/MW; % J/kg
g = 9.806; % m/s^2
gamma    = 1.25;
C_star = sqrt(R*T_c) / ( sqrt(gamma) * (2/(gamma+1))^((gamma+1)/(2*(gamma-1))) );

% Given: A_ratios = Ae/At array
A_ratios = linspace(30,100,1000);

% --- Solve for Mach on both branches (robust bracketing) ---
M_sup = zeros(size(A_ratios));
M_sub = zeros(size(A_ratios));

for i = 1:numel(A_ratios)
    A = A_ratios(i);

    f = @(M) (1./M) .* ...
        ( (2./(gamma+1)) .* (1 + (gamma-1)./2 .* M.^2) ) ...
        .^ ( (gamma+1)./(2*(gamma-1)) ) ...
        - A;
    M_sub(i) = fzero(f, [0.000001, 1]);
    M_sup(i) = fzero(f, [1, 10000]);    % supersonic root
end

% --- Pressure ratios from Mach (isentropic) ---
% p/p0 and p0/p for each branch
p_ratio_sup = (1 + (gamma-1)/2 .* M_sup.^2) .^ (-gamma/(gamma-1));
p_ratio_sub = (1 + (gamma-1)/2 .* M_sub.^2) .^ (-gamma/(gamma-1));

CF_ideal = sqrt( (2*gamma^2/(gamma-1)) .* (2/(gamma+1)).^((gamma+1)/(gamma-1)) ...
    .* (1 - (p_ratio_sup).^((gamma-1)/gamma)) ) ...
     + (p_ratio_sup) .* A_ratios;

u_e = sqrt((2*gamma/(gamma-1)*R*T_c).*(1-(p_ratio_sup).^((gamma-1)/gamma)));

T_e = T_c.*(p_ratio_sup).^((gamma-1)/gamma);

A_traget = 60;

CF_ideal_target = interp1(A_ratios, CF_ideal, A_traget);
u_e_target = interp1(A_ratios, u_e, A_traget);
T_e_target = interp1(A_ratios, T_e, A_traget);
p_ratio_target = interp1(A_ratios, p_ratio_sup, A_traget);
M_target = interp1(A_ratios, M_sup, A_traget);

I_sp_ideal = C_star*CF_ideal_target/g;

m_dot_ideal = F_ideal/(I_sp_ideal*g);

A_t = (m_dot_ideal*C_star)./Pc;

A_t_target = 4e-6;
Pc_target = interp1(A_t, Pc, A_t_target);
Pe = (2/3)*Pc_target*p_ratio_target;

V_c = m_dot_ideal .* t_r .* R .* T_c ./ Pc_target;
L_star = V_c/A_t_target;
L_c = 0.08; % want it to be 10 cm long
epsilon = L_star/L_c; %sizing factor also expansion ratio for subsonic region
A_c = A_t_target*epsilon;

% sanity check
F = m_dot_ideal*u_e_target + Pc_target*p_ratio_target*A_t_target*A_traget;

% Create table
Results = table( ...
    A_traget, ...
    Pc_target, ...
    Pe, ...
    T_e_target, ...
    'VariableNames', {'Expansion_Ratio', ...
                      'Chamber_Pressure_Pa', ...
                      'Exit_Pressure_Pa', ...
                      'Exit_Temperature_K'});


disp(Results)

%%
% matches? time for MoC to get a 2D contour of nozzle
rt = sqrt(A_t_target/pi());
N = 48;

[wall_pts, X, Y] = moc_minlen_net(gamma, M_target, rt, N);

fprintf('Wall Coordinates (x, y):\n');
for k = 1:size(wall_pts,1)
    fprintf('%.6f, %.6f\n', wall_pts(k,1), wall_pts(k,2));
end

% Plotting
figure('Color', 'w');
hold on; axis equal; grid on;

% Plot Characteristic Mesh (Red Lines)
% Fan Rays (C-): Connect (i,0) to (i,i)
for i = 1:N
    % Points for fixed i, varying j from 0 to i
    % indices in X/Y are (i, 1) to (i, i+1)
    plot(X(i, 1:i+1), Y(i, 1:i+1), 'r-', 'LineWidth', 0.5);
end

% Reflected Rays (C+): Connect (j,j) to (N,j) then to Wall
for j = 1:N
    % Kernel part: varying i from j to N for fixed j
    % indices in X/Y are (j:N, j+1)
    x_ray = X(j:N, j+1);
    y_ray = Y(j:N, j+1);

    % Extension to Wall
    % Ray j ends at wall index j+1 (since wall index 1 is throat start)
    x_ray(end+1) = wall_pts(j+1, 1);
    y_ray(end+1) = wall_pts(j+1, 2);

    plot(x_ray, y_ray, 'r-', 'LineWidth', 0.5);
end

% Plot Wall (Black Thick Line)
plot(wall_pts(:,1), wall_pts(:,2), 'k-', 'LineWidth', 2);

% Plot Axis (Black Line)
yline(0, 'k-', 'LineWidth', 1.5);

% Labels
xlabel('Axial Position x (m)');
ylabel('Radius r (m)');
title('Minimum Length Nozzle Design (MoC)');

% Set limits to show full nozzle
xlim([0, max(wall_pts(:,1))*1.05]);
ylim([0, max(wall_pts(:,2))*1.2]);

function [wall_pts, X, Y] = moc_minlen_net(gam, Me, Rt, N)
% gam  : gamma
% Me   : target exit Mach
% Rt   : throat radius (sets scale)
% N    : number of characteristic lines (30-80 typical)

% ---------- PM + inverse ----------
nuPM = @(M) sqrt((gam+1)/(gam-1))*atan( sqrt((gam-1)/(gam+1)*(M.^2-1)) ) ...
          - atan( sqrt(M.^2-1) );

nu_max = 0.5*pi*( sqrt((gam+1)/(gam-1)) - 1 );

invPM = @(nu_target) invPM_safe(nu_target, gam, nuPM, nu_max);

% muMach handles M=1 safely
    function val = muMach(M)
        if M <= 1
            val = pi/2;
        else
            arg = 1./M;
            if arg > 1
                arg = 1;
            end
            val = asin(arg);
        end
    end

% ---------- Target turning ----------
nu_e = nuPM(Me);
theta_max = 0.5*nu_e;                 % minimum-length (planar baseline)
dtheta = theta_max / N;

% Data Structure for Points
% We need random access by (i, j).
% Since i goes 1..N and j goes 0..i, we can use a Cell Array or simple structs
% Or just 2D arrays since variables are scalars.
% x(i, j+1) to handle 0-based j in 1-based MATLAB.
% Let's use 2D arrays.
% Rows i = 1..N. Cols j = 0..N.
% Index (i, j) in logic maps to (i, j+1) in MATLAB array.

sz = N + 1; % indices 1..N+1 for j=0..N
X   = zeros(N, sz);
Y   = zeros(N, sz);
Th  = zeros(N, sz);
Nu  = zeros(N, sz);
M   = zeros(N, sz);
Mu  = zeros(N, sz);

% Initialize Corner (i, j=0) -> MATLAB col 1
% State at corner for ray i: theta = i*dtheta, nu = theta.
for i = 1:N
    th_val = i * dtheta;
    nu_val = th_val;
    M_val  = invPM(nu_val);
    mu_val = muMach(M_val);

    X(i, 1)  = 0;
    Y(i, 1)  = Rt;
    Th(i, 1) = th_val;
    Nu(i, 1) = nu_val;
    M(i, 1)  = M_val;
    Mu(i, 1) = mu_val;
end

% Marching loop
% i (Fan Ray) from 1 to N
%   j (Refl Ray) from 1 to i
for i = 1:N
    for j = 1:i
        % MATLAB indices
        % Current Point P: (i, j) -> (i, j+1)
        % Left Point L: (i, j-1) -> (i, j)
        % Bottom Point B: (i-1, j) -> (i-1, j+1)

        idx_P = j + 1;
        idx_L = j;

        if j == i
            % Case 1: Axis Point (j == i)
            % Intersection of Fan Ray i (from L) with Axis.

            % L state
            x_L = X(i, idx_L);
            y_L = Y(i, idx_L);
            th_L = Th(i, idx_L);
            nu_L = Nu(i, idx_L);
            mu_L = Mu(i, idx_L);

            % Axis BC
            th_P = 0;
            y_P = 0;

            % Source term approximation (Planar: term=0)
            term_src = 0;

            K_minus_L = th_L - nu_L;
            % nu_P = - (K_minus_L + term_src)
            nu_P = -(K_minus_L + term_src);
            if nu_P < 0; nu_P = 0; end

            % Geometry
            m_minus = tan(th_L - mu_L);
            x_P = x_L + (y_P - y_L) / m_minus;

            M_P = invPM(nu_P);
            mu_P = muMach(M_P);

            % Store P
            X(i, idx_P) = x_P;
            Y(i, idx_P) = y_P;
            Th(i, idx_P) = th_P;
            Nu(i, idx_P) = nu_P;
            M(i, idx_P) = M_P;
            Mu(i, idx_P) = mu_P;

        else
            % Case 2: Interior Point (j < i)
            % L = (i, j-1) (From Fan Ray i) -> idx_L
            % B = (i-1, j) (From Refl Ray j) -> idx_B
            idx_B = j + 1;

            x_L = X(i, idx_L); y_L = Y(i, idx_L);
            th_L = Th(i, idx_L); nu_L = Nu(i, idx_L); mu_L = Mu(i, idx_L);

            x_B = X(i-1, idx_B); y_B = Y(i-1, idx_B);
            th_B = Th(i-1, idx_B); nu_B = Nu(i-1, idx_B); mu_B = Mu(i-1, idx_B);

            % Intersect lines
            m_plus = tan(th_B + mu_B);
            m_minus = tan(th_L - mu_L);

            x_P = (y_L - y_B + m_plus*x_B - m_minus*x_L) / (m_plus - m_minus);
            y_P = y_B + m_plus*(x_P - x_B);

            % Source terms (Planar: 0)
            RHS_plus = 0;
            RHS_minus = 0;

            J_plus = th_B + nu_B;
            J_minus = th_L - nu_L;

            sum_P = J_plus + RHS_plus;
            diff_P = J_minus + RHS_minus;

            th_P = 0.5 * (sum_P + diff_P);
            nu_P = 0.5 * (sum_P - diff_P);

            if nu_P < 0; nu_P = 0; end

            M_P = invPM(nu_P);
            mu_P = muMach(M_P);

            X(i, idx_P) = x_P;
            Y(i, idx_P) = y_P;
            Th(i, idx_P) = th_P;
            Nu(i, idx_P) = nu_P;
            M(i, idx_P) = M_P;
            Mu(i, idx_P) = mu_P;
        end
    end
end

% --- Wall Construction ---
% Trace from Point(N, j) (End of Refl Ray j in Kernel) to Wall
% Wall starts at (0, Rt).

% Wall points list.
wall_x = zeros(N+1, 1);
wall_y = zeros(N+1, 1);

% Start point
wall_x(1) = 0;
wall_y(1) = Rt;
current_th = theta_max;
% Actually we need to track current point.
wx = 0; wy = Rt; wth = theta_max;

for j = 1:N
    % Point P is on Fan Ray N, Refl Ray j.
    % Index in our arrays: P(N, j) -> row N, col j+1.
    idx_P = j + 1;

    x_P = X(N, idx_P);
    y_P = Y(N, idx_P);
    th_P = Th(N, idx_P);
    nu_P = Nu(N, idx_P);
    mu_P = Mu(N, idx_P);

    % Extension of Refl Ray j from P to Wall
    m_plus = tan(th_P + mu_P);

    % Wall Segment Target Theta
    % Linear decay
    th_wall_target = theta_max * (1.0 - j/N);

    % Wall Segment Slope (Average)
    th_avg = 0.5 * (wth + th_wall_target);
    m_wall = tan(th_avg);

    % Intersect
    x_W = (wy - y_P + m_plus*x_P - m_wall*wx) / (m_plus - m_wall);
    y_W = wy + m_wall*(x_W - wx);

    wall_x(j+1) = x_W;
    wall_y(j+1) = y_W;

    % Update current wall point
    wx = x_W;
    wy = y_W;
    wth = th_wall_target;
end

wall_pts = [wall_x, wall_y];

end

% -------- local helper --------
function M = invPM_safe(nu_target, gam, nuPM, nu_max)
    if ~isfinite(nu_target) || nu_target <= 1e-9
        M = 1.0;
        return;
    end
    if nu_target > nu_max
        nu_target = nu_max - 1e-9;
    end
    f = @(M) nuPM(M) - nu_target;

    % Robust fzero
    try
        M = fzero(f, [1, 1e5]);
    catch
        % Fallback if 1e5 is not enough or too much (unlikely for nu < nu_max)
        % Check sign at boundaries
        v1 = f(1);
        v2 = f(1e5);
        if v1*v2 > 0
             M = 1e5; % Saturation
        else
             M = fzero(f, 10); % Try single start point
        end
    end
end
