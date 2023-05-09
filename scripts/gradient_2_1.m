close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation data - second order system - three unknown parameters
% y/u = b0 / (s^2+a1s+a0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkg load control;
b0 = 2;
a1 = 2;
a0 = 1;
num = [0 0 b0];
den = [1 a1 a0];
ftx = tf(num,den);

N = 2000;
dt = 0.1;
T = 0:dt:N;
% input
%u = ones(1,N/dt+1);
u = 0.5*sin(T) + 1.2*cos(2*T);
%u = 0.1*cos((pi/7)*T) + 0.05*sin((pi/5)*T);
[Y,T] = lsim(ftx,u,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gradient method J instantaneous
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Erro(1)     = 0;
theta(:,1)  = [0;0;0];
gamma_1 = 1;
gamma_2 = 1;
gamma_3 = .5;
Gamma = [gamma_1 0 0 ;
         0 gamma_2 0;
         0 0 gamma_3];
alpha = 1;

%Filtros
f_z = tf([1 0 0], [1 2 1]);
z = lsim(f_z, Y, T);

f_phi_1 = tf([0 0 1],[1 2 1]);
f_phi_2 = tf([0 1 0], [1 2 1]);
f_phi_3 = f_phi_1;

phi_1 = lsim(f_phi_1,u,T);
phi_2 = lsim(f_phi_2,-Y,T);
phi_3 = lsim(f_phi_3,-Y,T);

for k = 1:N/dt
    Phi(:,k)   = [phi_1(k); phi_2(k); phi_3(k)];
    m2          = 1 + alpha*Phi(:,k)'*Phi(:,k);
    erro(k+1)   = (z(k) - theta(:,k)'*Phi(:,k))/m2;
    theta(:,k+1) = theta(:,k) + Gamma*dt*erro(k+1)*Phi(:,k);
end

subplot(3,1,1);
b0_plot = b0*ones(1, length(theta(1,:)));
plot(T, b0_plot, 'k-', 'LineWidth', 2, T, theta(1,:), 'k--', 'LineWidth', 2);
title('Gradient method - instantaneous cost function');
set(gca, "fontsize", 12);
axis([0, length(T)*dt, 0, 4]);
subplot(3,1,2);
a1_plot = a1*ones(1, length(theta(2,:)));
plot(T, a1_plot, 'k-', 'LineWidth', 2, T, theta(2,:), 'k--', 'LineWidth', 2);
set(gca, "fontsize", 12);
axis([0, length(T)*dt, 0, 4]);
subplot(3,1,3);
a0_plot = a0*ones(1, length(theta(3,:)));
plot(T, a0_plot, 'k-', 'LineWidth', 2, T, theta(3,:), 'k--', 'LineWidth', 2);
set(gca, "fontsize", 12);
axis([0, length(T)*dt, 0, 4]);
legend('nominal value', 'Parameter being identified: a', "location", "southeast");
xlabel('Samples');
%ylabel('Estimating a...');
%subplot(2,1,2);
%plot(T, erro, 'k-', 'LineWidth', 2);
%set(gca, "fontsize", 12);
%xlabel('Samples');
%ylabel('Estimation error');
