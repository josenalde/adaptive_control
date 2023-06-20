close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation data - first order system - two unknown parameters
% y/u = b0 / (s + a1s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkg load control;
kp = 0.5;
a1 = 1;
km = 2;
am1 = 2;

dt = 0.1;
tMax = 100;

gamma_1 = 0.5;
gamma_2 = 1;

t = 0;


u = zeros(1, tMax/dt + 2);
tplot = zeros(1, tMax/dt + 2);
y = zeros(1, tMax/dt + 2);
ym = zeros(1, tMax/dt + 2);
k = 2;
e = zeros(1, tMax/dt + 2);
theta_1 = zeros(1, tMax/dt + 2);
theta_2 = zeros(1, tMax/dt + 2);
theta_1_star = ((a1 - am1) / kp) *ones(1, tMax/dt + 2);
theta_2_star = (km/kp)*ones(1, tMax/dt + 2);

while (t < tMax)

   y(k) = y(k-1) + dt*(-a1*y(k-1) + kp*u(k-1));
   r(k) = 1;
   %r(k) = 1 + 0.3*sin(t) + 0.3*cos(2*t);
   ym(k) = ym(k-1) + dt*(-am1*ym(k-1) + km*r(k));
   %u(k) = theta_1_star*y(k) + theta_2_star*r(k);
   e(k) = y(k) - ym(k);
   % leis adaptativas Lyapunov
   theta_1(k) = theta_1(k-1) + dt*(-gamma_1 * e(k) * y(k));
   theta_2(k) = theta_2(k-1) + dt*(-gamma_2 * e(k) * r(k));
   u(k) = theta_1(k)*y(k) + theta_2(k)*r(k);
   t = t + dt;
   tplot(k) = t;
   k = k + 1;
endwhile

plot(tplot(1:length(tplot)-2), y(1:length(tplot)-2), 'k', 'LineWidth', 2, tplot(1:length(tplot)-2), ym(1:length(tplot)-2), 'b', 'LineWidth', 2);
figure;
plot(tplot(1:length(tplot)-2), theta_1(1:length(tplot)-2), 'k', 'LineWidth', 2, tplot(1:length(tplot)-2), theta_1_star(1:length(tplot)-2), 'b--', 'LineWidth', 2);
figure;
plot(tplot(1:length(tplot)-2), theta_2(1:length(tplot)-2), 'k', 'LineWidth', 2, tplot(1:length(tplot)-2), theta_2_star(1:length(tplot)-2), 'b--', 'LineWidth', 2);
clear;

