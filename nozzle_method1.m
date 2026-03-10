% Nozzle Contour Generation - Method 1: Rao's Parabolic Approximation
%
% Summary of Logic:
% This script generates a 2D nozzle contour using the standard Rao Parabolic Approximation method.
% 1. It calculates the necessary nozzle parameters (Area Ratio, Throat Radius Rt, Exit Radius Re)
%    based on the given inputs (Thrust, Pressure, Temperature, Gamma).
% 2. It defines the nozzle geometry in three sections:
%    - Convergent Section: A circular arc leading to the throat.
%    - Throat Section: A circular arc downstream of the throat (radius 0.382 * Rt).
%    - Divergent Bell Section: A parabola fitting between the throat section end and the exit radius.
%      The parabola is defined by the initial wall angle (theta_n) and the exit wall angle (theta_e).
% 3. The script solves for the parabolic coefficients and generates coordinate points.
% 4. Finally, it plots the full 2D contour and outputs the coordinates to the command window.

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

% Design Point Selection (Using target values from original code)
A_ratio_target = 60;
A_t_target = 4e-6; % m^2

% Calculate Exit Radius
Rt = sqrt(A_t_target/pi);
Re = sqrt(A_ratio_target) * Rt;

% --- Rao Nozzle Geometry Parameters ---
% Standard parabolic bell nozzle parameters (approximate optimum)
% Length of nozzle: often defined as percentage of a 15-degree conical nozzle
% L_cone_15 = (Re - Rt) / tan(deg2rad(15));
% L_nozzle = 0.8 * L_cone_15; % 80% bell nozzle

theta_n = 28 * (pi/180); % Initial parabola angle (inflection point angle), typical 20-30 deg
theta_e = 6 * (pi/180);  % Exit angle, typical 2-8 deg

% Coordinate System: Throat is at (0, Rt)
% Throat Downstream Arc: Radius Rw = 0.382 * Rt (Standard)
Rw = 0.382 * Rt;

% End of Throat Arc (Inflection Point N)
x_n = Rw * sin(theta_n);
y_n = Rt + Rw * (1 - cos(theta_n));

% Length of Bell (L)
% We can determine L based on fitting a parabola y = ax^2 + bx + c
% Constraints:
% 1. At x = x_n, y = y_n
% 2. At x = x_n, dy/dx = tan(theta_n)
% 3. At x = L, y = Re
% 4. At x = L, dy/dx = tan(theta_e)
% This is an over-constrained system for a simple quadratic parabola if L is fixed.
% Standard Rao method actually constructs the parabola based on these boundary conditions to FIND L.
% Parabola: x(y) or y(x). Usually y(x) = a*x^2 + b*x + c is used.
% Let's use coordinate transformation so x_n is 0 for the parabola segment.
% y = a*x'^2 + b*x' + c  where x' = x - x_n.
% At x'=0 (Inflection): y = y_n, dy/dx' = tan(theta_n)
%   -> c = y_n
%   -> b = tan(theta_n)
% At x'=Ln (Exit): y = Re, dy/dx' = tan(theta_e)
%   -> Re = a*Ln^2 + b*Ln + c
%   -> tan(theta_e) = 2*a*Ln + b
% Solve for a and Ln.
% From (2): a = (tan(theta_e) - b) / (2*Ln)
% Subst into (1): Re = [(tan(theta_e) - b)/(2*Ln)] * Ln^2 + b*Ln + c
% Re = (tan(theta_e) - b)*Ln/2 + b*Ln + c
% Re - c = Ln/2 * (tan(theta_e) - b + 2*b)
% Re - y_n = Ln/2 * (tan(theta_e) + tan(theta_n))
% Ln = 2 * (Re - y_n) / (tan(theta_e) + tan(theta_n))
% a = (tan(theta_e) - tan(theta_n)) / (2*Ln)

Ln = 2 * (Re - y_n) / (tan(theta_e) + tan(theta_n));
a  = (tan(theta_e) - tan(theta_n)) / (2*Ln);
b  = tan(theta_n);
c  = y_n;

% --- Generate Points ---
N_points = 100;

% 1. Throat Arc Section
theta_arc = linspace(0, theta_n, round(N_points/3));
x_arc = Rw * sin(theta_arc);
y_arc = Rt + Rw * (1 - cos(theta_arc));

% 2. Bell Parabola Section
x_par_local = linspace(0, Ln, round(2*N_points/3));
x_bell = x_n + x_par_local;
y_bell = a * x_par_local.^2 + b * x_par_local + c;

% Combine
x_nozzle = [x_arc, x_bell];
y_nozzle = [y_arc, y_bell];

% --- Output ---
fprintf('Method 1: Rao Parabolic Approximation\n');
fprintf('Expansion Ratio: %.2f\n', A_ratio_target);
fprintf('Nozzle Length:   %.4f m\n', x_nozzle(end));
fprintf('Exit Radius:     %.4f m\n', y_nozzle(end));
fprintf('Wall Coordinates (x, y):\n');
for i = 1:10:length(x_nozzle) % Print sparse sample
    fprintf('%.6f, %.6f\n', x_nozzle(i), y_nozzle(i));
end

% --- Plotting ---
figure('Name', 'Method 1: Rao Parabolic');
hold on; axis equal; grid on;

% Fill Interior
fill([x_nozzle, fliplr(x_nozzle)], [y_nozzle, -fliplr(y_nozzle)], [0.9 0.9 1.0], 'EdgeColor', 'none');

% Plot Upper and Lower Walls
plot(x_nozzle, y_nozzle, 'k-', 'LineWidth', 2.5);
plot(x_nozzle, -y_nozzle, 'k-', 'LineWidth', 2.5);

% Plot Axis
yline(0, 'k-.', 'LineWidth', 1.5);
xline(0, 'k--', 'LineWidth', 1);

% --- Annotations ---

% 1. Radii (Rt and Re)
plot([0 0], [0 Rt], 'r-', 'LineWidth', 1.5); % Rt Line
plot([x_nozzle(end) x_nozzle(end)], [0 Re], 'r-', 'LineWidth', 1.5); % Re Line

% Rt Label: Left and Down
text(0, Rt/3, sprintf('R_t = %.4f m', Rt), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', 90, 'FontSize', 10, 'Color', 'r');
% Re Label: More Right
text(x_nozzle(end), Re/3, sprintf('R_e = %.4f m', Re), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Rotation', 90, 'FontSize', 10, 'Color', 'r');

% 2. Inflection Point (xn, yn)
plot(x_n, y_n, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
text(x_n, y_n*1.1, sprintf('Inflection\n\\theta_n = %.1f^o', rad2deg(theta_n)), 'HorizontalAlignment', 'center', 'FontSize', 10);

% 3. Exit Point Angle
plot(x_nozzle(end), y_nozzle(end), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
text(x_nozzle(end), y_nozzle(end)*1.1, sprintf('Exit\n\\theta_e = %.1f^o', rad2deg(theta_e)), 'HorizontalAlignment', 'center', 'FontSize', 10);

% 4. Throat Label
text(-0.001, Rt*1.2, 'Throat Plane', 'HorizontalAlignment', 'right', 'FontSize', 10);

% Axes Labels and Title
xlabel('Axial Position x (m)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Radius r (m)', 'FontSize', 12, 'FontWeight', 'bold');
title({'2D Nozzle Contour - Rao Parabolic Approximation'; sprintf('Area Ratio = %.1f', A_ratio_target)}, 'FontSize', 14);

% Limits
xlim([-0.005, x_nozzle(end)*1.1]);
ylim([-Re*1.3, Re*1.3]);
