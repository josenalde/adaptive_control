%script adaptado de https://edisciplinas.usp.br/course/view.php?id=81578#section-1 #

close all
clear all

pkg load control;
% parâmetros nominais (ou seja, do sistema real) de onde se medirá IO
% sistema
m  = 15;
b = 0.2;
km = 2;
num = [0 0 1];
den = [m b (km-1)];
ftx  = tf(num,den);

% dados para simulação
N = 150;
dt = 0.1;
T = 0:dt:N;
u = ones(1,N/dt+1);
u = 10*sin(2*T) + 7*cos(3.6*T);
%Simula o sistema para a entrada u degrau unitario (saída do sistema nominal - "real")
[Y,T] = lsim(ftx,u,T);

%plot(T,Y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% algoritmo do gradiente para identificar parametros
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parâmetros
theta(:,1)  = [0.5; 0.3]; % beta, k-1
rho(1) = 1/12;
gamma_1 = 20;
gamma_2 = 5;
gamma_rho = 0.5;
Gamma = [gamma_1 0; 0 gamma_2];
alpha = 1; % para o sinal de normalização ns

% z = dot(y) / Lambda
% Componentes de z com Lambda = s^2 + 2s + 1
f_z = tf([1 0 0],[1 2 1]); % s^2 / s^2+2s+1
f_phi_1 = tf([0 1 0],[1 2 1]); % s / s^2+2s+1
f_phi_2 = tf([0 0 1],[1 2 1]); % s / s^2+2s+1
f_u = tf([0 0 1],[1 2 1]);

z = lsim(f_z, Y,T);
phi_1 = lsim(f_phi_1,-Y,T);
phi_2 = lsim(f_phi_2,-Y,T);
z0 = lsim(f_u, u,T);

for k = 1:N/dt
    Phi(:,k)   = [phi_1(k); phi_2(k)];
    m2          = 1 + alpha*Phi(:,k)'*Phi(:,k);
    zeta(k+1) = theta(:,k)'*Phi(:,k) + z0(k);
    erro(k+1)   = (z(k) - rho(k) .* zeta(k+1))/m2;
    rho(k+1) = rho(k) + gamma_rho*erro(k+1)*zeta(k+1);
    theta(:,k+1) = theta(:,k) + Gamma*dt*erro(k+1)*Phi(:,k);
    % Projecao
##    if ((theta(1,k) > 1) || ((theta(1,k) == 1) && (gamma_1 * erro(k+1) * phi_1 >= 0)))
##      disp('entrei no if b');
##      theta(1,k+1) = theta(1,k) + gamma_1 * dt * erro(k+1) * phi_1(k);
##    elseif ((theta(1,k) == 1 && (gamma_1 * erro(k+1) * phi_1(k) < 0)))
##      theta(1,k+1) = theta(1,k);
##      disp('entrei no else b');
##    endif
##
##    if ((theta(2,k) > -2 && theta(2,k) < 20) || ((theta(2,k) == -2) && (gamma_2 * erro(k+1) * phi_2(k) >= 0)) || ((theta(2,k) == 20) && (gamma_2 * erro(k+1) * phi_2(k) <= 0)))
##      disp('entrei no if a');
##      theta(2,k+1) = theta(2,k) + gamma_2 * dt * erro(k+1) * phi_2(k);
##    elseif ((theta(2,k) == -2 && (gamma_2 * erro(k+1) * phi_2(k) < 0)) || (theta(2,k) == 20 && (gamma_2 * erro(k+1) * phi_2(k) > 0)))
##      theta(2,k+1) = theta(2,k);
##      disp('entrei no else a');
##    endif
end

subplot(3,1,1);
m_plot = m*ones(1, length(theta(1,:)));
%ajustando
mt = 1 ./ rho;
plot(T, m_plot, 'k-', 'LineWidth', 2, T, mt, 'k--', 'LineWidth', 2);
title('Gradient method - projection');
set(gca, "fontsize", 12);
%axis([0, length(T)*dt, 0, 4]);
subplot(3,1,2);
b_plot = b*ones(1, length(theta(1,:)));
plot(T, b_plot, 'k-', 'LineWidth', 2, T, theta(1,:), 'k--', 'LineWidth', 2);
set(gca, "fontsize", 12);
%axis([0, length(T)*dt, 0, 4]);
subplot(3,1,3);
k_plot = km*ones(1, length(theta(2,:)));
kt = theta(2,:)  + 1;
plot(T, k_plot, 'k-', 'LineWidth', 2, T, kt, 'k--', 'LineWidth', 2);
set(gca, "fontsize", 12);

xlabel('Samples');

