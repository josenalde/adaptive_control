%script adaptado de https://edisciplinas.usp.br/course/view.php?id=81578#section-1 #

close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%
% dados de simulacao
% \dot{x}(t) = -ax(t)+4u(t)
% x = (b/(s+a))u, com b conhecido
%%%%%%%%%%%%%%%%%%%%%%%%%%
pkg load control;
% parâmetros nominais (ou seja, do sistema real) de onde se medirá IO
% sistema
a  = 3;
b = 4;
num = [0 b];
den = [1 a];
ftx  = tf(num,den);

% dados para simulação
N = 250;
dt = 0.1;
T = 0:dt:N;
u = ones(1,N/dt+1);

%Simula o sistema para a entrada u degrau unitario (saída do sistema nominal - "real")
[Y,T] = lsim(ftx,u,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% algoritmo do gradiente para identificar parametros
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parâmetros
theta(1)     = 0; %valor inicial de theta = a
gamma        = 0.1; %ganho adaptativo
alpha = 1; % para o sinal de normalização ns

% z = dot(y) / Lambda - (b /Lambda) u
% Componentes de z com Lambda = s + 2
f_dot_y = tf([1 0],[1 2]); % s / s+2
f_phi_1 = tf([0 1],[1 2]); % 1 / s+2
f_u = tf(4,[1 2]);

Y2 = lsim(f_dot_y,Y,T);
Y3 = lsim(f_u,u,T);
z = Y2 - Y3;

Y1 = lsim(f_phi_1,Y,T);

for k = 1:N/dt
    phi_1(k)   = -Y1(k);
    m2         = 1 + alpha*(phi_1(k)'*phi_1(k));
    erro(k+1) = (z(k) - theta(k)'*phi_1(k)) / m2;
    theta(k+1) = theta(k) + gamma*dt*erro(k+1)*phi_1(k);
end

subplot(2,1,1);
a_plot = a*ones(1, length(theta));
plot(T, a_plot, 'k-', 'LineWidth', 2, T, theta, 'k--', 'LineWidth', 2);
set(gca, "fontsize", 12);
axis([0, length(T)*dt, 0, 4]);
title('Gradient method - instantaneous cost function');
legend('nominal value', 'Parameter being identified: a', "location", "southeast");
xlabel('Samples');
ylabel('Estimating a...');
subplot(2,1,2);
plot(T, erro, 'k-', 'LineWidth', 2);
set(gca, "fontsize", 12);
xlabel('Samples');
ylabel('Estimation error');

