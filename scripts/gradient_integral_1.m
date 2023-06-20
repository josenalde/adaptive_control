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

dt = 0.1;
tMax = 100;

N = tMax/dt;
T = 0:dt:N;
% input
%u = ones(1,N/dt+1);

u = 0.5*sin(T) + 1.2*cos(2*T);
%u = 0.1*cos((pi/7)*T) + 0.05*sin((pi/5)*T);
[Y,T] = lsim(ftx,u,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gradient method J integral
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Erro(1)     = 0;
theta(:,1)  = [0;0];
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

phi_1 = lsim(f_phi_1, u, T);
phi_2 = lsim(f_phi_2, -Y, T);

k = 2;

theta_1_star = 0.5;
theta_2_star = 1;

theta_1_star = theta_1_star * ones(1, length(T));
theta_2_star = theta_2_star * ones(1, length(T));

beta = 0.1; %para função de custo integral

q(:,1) = [0;0];
R(:,:,1) = [0 0; 0 0];

for k = 1:N/dt
   Phi(:,k)   = [phi_1(k); phi_2(k)];
   m2          = 1 + alpha*Phi(:,k)'*Phi(:,k);
   q(:,k+1) = q(:,k) + dt*(-beta * q(:,k) - ((z(k) * Phi(:,k)) / m2));
   R(:,:,k+1) = R(:,:,k) + dt*(-beta * R(:,:,k) + ((Phi(:,k)*Phi(:,k)')/m2));
   theta(:,k+1) = theta(:,k) + dt*(-Gamma*(R(:,:,k+1)*theta(:,k) + q(:,k+1)));
endfor

plot(T, theta(1,:), 'k', 'LineWidth', 2, T, theta_1_star, 'b--', 'LineWidth', 2);
figure;
plot(T, theta(2,:), 'k', 'LineWidth', 2, T, theta_2_star, 'b--', 'LineWidth', 2);

clear;

