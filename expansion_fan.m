% Summary:
% This script computes and plots the velocity and pressure distributions
% for a Prandtl-Meyer expansion fan flow using the Method of Characteristics.
% It solves for the flow properties (velocity u and pressure p) as a function of
% the similarity variable eta = x/(a1*t).

gamma = 1.4;
U = 0.3;                % example u2/a1
eta = linspace(-2, 1, 2000);

eta_head = -1;
eta_tail = -1 + (gamma+1)/2 * U;

u = zeros(size(eta));
p = ones(size(eta));

% fan indices
I = (eta >= eta_head) & (eta <= eta_tail);

u(I) = (2/(gamma+1)) * (1 + eta(I));     % u/a1 in fan
u(eta > eta_tail) = U;                   % region 2

a_over_a1 = ones(size(eta));
a_over_a1(I) = (2 - (gamma-1)*eta(I))/(gamma+1);
a_over_a1(eta > eta_tail) = 1 - (gamma-1)/2 * U;

p(I) = a_over_a1(I).^(2*gamma/(gamma-1));
p(eta > eta_tail) = a_over_a1(eta > eta_tail).^(2*gamma/(gamma-1));

% Plotting Velocity Distribution
figure;
plot(eta, u, 'LineWidth', 2);
grid on;
title('Velocity Distribution in Expansion Fan', 'FontSize', 16);
xlabel('x/(a_1 t)', 'FontSize', 14);
ylabel('u/a_1', 'FontSize', 14);
set(gca, 'FontSize', 12);

% Plotting Pressure Distribution
figure;
plot(eta, p, 'LineWidth', 2);
grid on;
title('Pressure Distribution in Expansion Fan', 'FontSize', 16);
xlabel('x/(a_1 t)', 'FontSize', 14);
ylabel('p/p_1', 'FontSize', 14);
set(gca, 'FontSize', 12);
