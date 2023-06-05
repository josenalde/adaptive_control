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
a  = 10;
b = 5;
num = [0 b];
den = [1 a];
ftx  = tf(num,den);

% dados para simulação
N = 350;
dt = 0.1;
T = 0:dt:N;
%u = ones(1,N/dt+1);
u = 5*sin(2*T);
%Simula o sistema para a entrada u degrau unitario (saída do sistema nominal - "real")
[Y,T] = lsim(ftx,u,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% algoritmo do gradiente para identificar parametros
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parâmetros
theta(:,1)  = [2; 2];
% theta_1 = b >= 1
% theta_2 = a <= 20 && a >= -2
gamma_1 = 1.4;
gamma_2 = 4;
Gamma = [gamma_1 0; 0 gamma_2];
alpha = 1; % para o sinal de normalização ns

% z = dot(y) / Lambda
% Componentes de z com Lambda = s + 2
f_dot_y = tf([1 0],[1 2]); % s / s+2
f_phi_1 = tf([0 1],[1 2]); % 1 / s+2
f_phi_2 = tf([0 1],[1 2]); % 1 / s+2

z = lsim(f_dot_y,Y,T);
phi_1 = lsim(f_phi_1,u,T);
phi_2 = lsim(f_phi_2,-Y,T);

phi = [phi_1; phi_2];

for k = 1:N/dt
    Phi(:,k)   = [phi_1(k); phi_2(k)];
    m2          = 1 + alpha*Phi(:,k)'*Phi(:,k);
    erro(k+1)   = (z(k) - theta(:,k)'*Phi(:,k))/m2;
    %m2_b = 1 + alpha*phi_1(k)*phi_1(k);
    %m2_a = 1 + alpha*phi_2(k)*phi_2(k);
    %erro_b(k+1) = (z(k) - theta(1,k)'*phi_1(k))/m2_b;
    %erro_a(k+1) = (z(k) - theta(2,k)'*phi_2(k))/m2_a;
    %theta(:,k+1) = theta(:,k) + Gamma*dt*erro(k+1)*Phi(:,k);
    % Projecao
    if ((theta(1,k) > 1) || ((theta(1,k) == 1) && (gamma_1 * erro(k+1) * phi_1 >= 0)))
      theta(1,k+1) = theta(1,k) + gamma_1 * dt * erro(k+1) * phi_1(k);
    elseif ((theta(1,k) == 1 && (gamma_1 * erro(k+1) * phi_1(k) < 0)))
      theta(1,k+1) = theta(1,k);
    endif

    if ((theta(2,k) > -2 && theta(2,k) < 20) || ((theta(2,k) == -2) && (gamma_2 * erro(k+1) * phi_2(k) >= 0)) || ((theta(2,k) == 20) && (gamma_2 * erro(k+1) * phi_2(k) <= 0)))
      theta(2,k+1) = theta(2,k) + gamma_2 * dt * erro(k+1) * phi_2(k);
    elseif ((theta(2,k) == -2 && (gamma_2 * erro(k+1) * phi_2(k) < 0)) || (theta(2,k) == 20 && (gamma_2 * erro(k+1) * phi_2(k) > 0)))
      theta(2,k+1) = theta(2,k);
    endif
end

subplot(2,1,1);
b_plot = b*ones(1, length(theta(1,:)));
plot(T, b_plot, 'k-', 'LineWidth', 2, T, theta(1,:), 'k--', 'LineWidth', 2);
title('Gradient method - projection');
set(gca, "fontsize", 12);
%axis([0, length(T)*dt, 0, 4]);
subplot(2,1,2);
a_plot = a*ones(1, length(theta(2,:)));
plot(T, a_plot, 'k-', 'LineWidth', 2, T, theta(2,:), 'k--', 'LineWidth', 2);
set(gca, "fontsize", 12);
%axis([0, length(T)*dt, 0, 4]);
xlabel('Samples');

