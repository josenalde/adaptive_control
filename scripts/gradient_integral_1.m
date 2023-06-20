close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation data - first order system - two unknown parameters
% y/u = b0 / (s + a1s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkg load control;
kp = 0.5;
a1 = 1;
num = [0 kp];
den = [1 a1];
ftx = tf(num,den);
km = 2;
am1 = 2;

dt = 0.1;
tMax = 500;

N = tMax/dt;
T = 0:dt:N;
% input
%u = ones(1,N/dt+1);
%r = ones(1,N/dt+1);
r = 1 + 0.3*sin(T) + 0.3*cos(2*T);
%r = 1 + (0.5*sin(T) + 1.2*cos(2*T));
%u = 0.1*cos((pi/7)*T) + 0.05*sin((pi/5)*T);
[Y,T] = lsim(ftx,r,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gradient method J instantaneous
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Erro(1)     = 0;
theta(:,1)  = [-1.0;3];
gamma_1 = 0.5;
gamma_2 = 0.1;
Gamma = [gamma_1 0;
         0 gamma_2];
alpha = 1;

%Filtros
f_z = tf([1 0], [1 2]);
z = lsim(f_z, Y, T);

f_phi_1 = tf([0 1],[1 2]);
f_phi_2 = tf([0 1], [1 2]);

%phi_1 = lsim(f_phi_1, r, T);
phi_2 = lsim(f_phi_2, -Y, T);

t = 0;


u = zeros(1, tMax/dt + 2);
tplot = zeros(1, tMax/dt + 2);
y = zeros(1, tMax/dt + 2);
ym = zeros(1, tMax/dt + 2);
k = 2;
e = zeros(1, tMax/dt + 2);
theta_1 = zeros(1, tMax/dt + 2);
theta_2 = zeros(1, tMax/dt + 2);
phi_1(1) = 0;


theta_1_star = ((a1 - am1) / kp) *ones(1, tMax/dt + 2);
theta_2_star = (km/kp)*ones(1, tMax/dt + 2);

beta = 0.1; %para função de custo integral

q(:,1) = [0;0];
R(:,:,1) = [0 0; 0 0];


while (t < tMax)

   y(k) = y(k-1) + dt*(-a1*y(k-1) + kp*u(k-1));
   r(k) = 1;
   %r(k) = 0.1*cos((pi/7)*t) + 0.05*sin((pi/5)*t);
   %r(k) = 1 + 0.3*sin(t) + 0.3*cos(2*t);
   ym(k) = ym(k-1) + dt*(-am1*ym(k-1) + km*r(k));
   %u(k) = theta_1_star*y(k) + theta_2_star*r(k);
   e(k) = y(k) - ym(k);
   % leis adaptativas Lyapunov
   theta_1(k) = theta_1(k-1) + dt*(-gamma_1 * e(k) * y(k));
   theta_2(k) = theta_2(k-1) + dt*(-gamma_2 * e(k) * r(k));
   % leis adaptativas gradiente integral
   phi_1(k) = phi_1(k-1) + dt*(-2*phi_1(k-1) + u(k-1));
   Phi(:,k)   = [phi_1(k); phi_2(k)];
   m2          = 1 + alpha*Phi(:,k)'*Phi(:,k);
   q(:,k) = q(:,k-1) + dt*(-beta * q(:,k-1) - ((z(k) * Phi(:,k)) / m2));
   R(:,:,k) = R(:,:,k-1) + dt*(-beta * R(:,:,k-1) + ((Phi(:,k)*Phi(:,k)')/m2));
   theta(:,k) = theta(:,k-1) + dt*(-Gamma*(R(:,:,k)*theta(:,k-1) + q(:,k)));
   %u(k) = theta_1(k)*y(k) + theta_2(k)*r(k);
   u(k) = theta(1,k)*y(k) + theta(2,k)*r(k);
   t = t + dt;
   tplot(k) = t;
   k = k + 1;
endwhile

plot(tplot(1:length(tplot)-2), y(1:length(tplot)-2), 'k', 'LineWidth', 2, tplot(1:length(tplot)-2), ym(1:length(tplot)-2), 'b', 'LineWidth', 2);
figure;
%plot(tplot(1:length(tplot)-2), theta_1(1:length(tplot)-2), 'k', 'LineWidth', 2, tplot(1:length(tplot)-2), theta_1_star(1:length(tplot)-2), 'b--', 'LineWidth', 2);
plot(tplot, theta(1,:), 'k', 'LineWidth', 2, tplot(1:length(tplot)-2), theta_1_star(1:length(tplot)-2), 'b--', 'LineWidth', 2);
figure;
%plot(tplot(1:length(tplot)-2), theta_2(1:length(tplot)-2), 'k', 'LineWidth', 2, tplot(1:length(tplot)-2), theta_2_star(1:length(tplot)-2), 'b--', 'LineWidth', 2);
plot(tplot, theta(2,:), 'k', 'LineWidth', 2, tplot(1:length(tplot)-2), theta_2_star(1:length(tplot)-2), 'b--', 'LineWidth', 2);
##for k = 1:N/dt
##    omega(:,k) = [Y(k); r(k)];
##    Phi(:,k)   = [phi_1(k); phi_2(k)];
##    m2          = 1 + alpha*Phi(:,k)'*Phi(:,k);
##    erro(k+1)   = (z(k) - theta(:,k)'*Phi(:,k))/m2;
##    theta(:,k+1) = theta(:,k) + Gamma*dt*erro(k+1)*Phi(:,k);
##end
##
##subplot(2,1,1);
##b0_plot = b0*ones(1, length(theta(1,:)));
##plot(T, b0_plot, 'k-', 'LineWidth', 2, T, theta(1,:), 'k--', 'LineWidth', 2);
##title('Gradient method - integral cost function');
##set(gca, "fontsize", 12);
##axis([0, length(T)*dt, 0, 1]);
##subplot(2,1,2);
##a1_plot = a1*ones(1, length(theta(2,:)));
##plot(T, a1_plot, 'k-', 'LineWidth', 2, T, theta(2,:), 'k--', 'LineWidth', 2);
##set(gca, "fontsize", 12);
##axis([0, length(T)*dt, 0, 1.5]);
##xlabel('Samples');
clear;

