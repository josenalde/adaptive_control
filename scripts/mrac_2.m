%by Prof. Josenalde Oliveira, PEM - UFRN @2023
% MRAC order 2

close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation data - second order system
% y/u = b1 / (s^2 + a1s + a2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkg load control;
pkg load signal;

% system

kp = 1;
b1 = 1;
a1 = -2;
a2 = 1;

gs = tf([0 kp b1], [1 a1 a2]);

% reference model

km = 1;
bm1 = 2;
am1 = 4;
am2 = 3;

ms = tf([0 km bm1], [1 am1 am2]);

[A, B, C, D] = tf2ss(gs);
[Am, Bm, Cm, Dm] = tf2ss(ms);

dt = 0.1;
tMax = 50;

lambda1 = 2; %for IO filters v1,v2
gamma = 2; %for g vector used for IO filter's calculation

gamma1 = 1; % adaptative gain for \theta_v11
gamma2= 2; % adaptive gain for \theta_2
gamma3 = 1; % adaptative gain for \theta_v21
gamma4 = 10; % adaptive gain for \theta_4

% initializations
t = 0;
u = zeros(1, tMax/dt + 2);
r = ones(1, tMax/dt + 2);

y = zeros(1, tMax/dt + 2);
ym = zeros(1, tMax/dt + 2);
e = zeros(1, tMax/dt + 2);
theta_v1 = zeros(1, tMax/dt + 2);
theta_2 = zeros(1, tMax/dt + 2);
theta_v2 = zeros(1, tMax/dt + 2);
theta_4 = zeros(1, tMax/dt + 2);

theta_v1_star = ((bm1 - b1) / gamma)*ones(1, tMax/dt + 2);
theta_2_star = ((a1 - am1) / kp) *ones(1, tMax/dt + 2);
theta_v2_star = ((a2 - am2 - kp*((a1 - am1) / kp)*bm1)/(kp*gamma))*ones(1, tMax/dt + 2);
theta_4_star = (km/kp)*ones(1, tMax/dt + 2);

tplot = zeros(1, tMax/dt + 2);
k = 2;

X(:,1) = [0;0];
Xm(:,1) = [0;0];
v1 = zeros(1, tMax/dt + 2);
v2 = zeros(1, tMax/dt + 2);

g = gamma;

while (t < tMax)

   dotX = A*X(:,k-1) + B*u(k-1);
   X(:,k) = X(:,k-1) + dt*(dotX);
   y(:,k) = C*X(:,k);

   dotXm = Am*Xm(:,k-1) + Bm*r(k-1);
   Xm(:,k) = Xm(:,k-1) + dt*(dotXm);
   ym(:,k) = Cm*Xm(:,k);

   dotv1 = -lambda1*v1(k-1) + g*u(k-1);
   v1(k) = v1(k-1) + dt*(dotv1);

   dotv2 = -lambda1*v2(k-1) + g*y(k);
   v2(k) = v2(k-1) + dt*(dotv2);

   e(k) = y(k) - ym(k);

   % Lyapunov-based laws
   theta_v1(k) = theta_v1(k-1) + dt*(-gamma1 * e(k) * v1(k));
   theta_2(k) = theta_2(k-1) + dt*(-gamma2 * e(k) * y(k));
   theta_v2(k) = theta_v2(k-1) + dt*(-gamma3 * e(k) * v2(k));
   theta_4(k) = theta_4(k-1) + dt*(-gamma4 * e(k) * r(k));

   u(k) = theta_v1(k)*v1(k) + theta_2(k)*y(k) + theta_v2(k)*v2(k) + theta_4(k)*r(k);
   t = t + dt;
   tplot(k) = t;
   k = k + 1;
endwhile
##
subplot(2,1,1);
plot(tplot(1:length(tplot)-2), y(1,1:length(tplot)-2), 'k', 'LineWidth', 2, tplot(1:length(tplot)-2), ym(1,1:length(tplot)-2), 'k', 'LineWidth', 2);
title('Plant output');
subplot(2,1,2);
plot(tplot(1:length(tplot)-2), u(1,1:length(tplot)-2), 'k', 'LineWidth', 2);
title('Control signal');
xlabel('Time (s)');
figure;
subplot(4,1,1);
plot(tplot(1:length(tplot)-2), theta_v1(1:length(tplot)-2), 'k', 'LineWidth', 2, tplot(1:length(tplot)-2), theta_v1_star(1:length(tplot)-2), 'b--', 'LineWidth', 2);
subplot(4,1,2);
plot(tplot(1:length(tplot)-2), theta_2(1:length(tplot)-2), 'k', 'LineWidth', 2, tplot(1:length(tplot)-2), theta_2_star(1:length(tplot)-2), 'b--', 'LineWidth', 2);
subplot(4,1,3);
plot(tplot(1:length(tplot)-2), theta_v2(1:length(tplot)-2), 'k', 'LineWidth', 2, tplot(1:length(tplot)-2), theta_v2_star(1:length(tplot)-2), 'b--', 'LineWidth', 2);
subplot(4,1,4);
plot(tplot(1:length(tplot)-2), theta_4(1:length(tplot)-2), 'k', 'LineWidth', 2, tplot(1:length(tplot)-2), theta_4_star(1:length(tplot)-2), 'b--', 'LineWidth', 2);


