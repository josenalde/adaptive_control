%script adaptado de https://edisciplinas.usp.br/course/view.php?id=81578#section-1 #

close all
clear all

%para uso no octave necessário linha abaixo
pkg load control;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dados de simulacao
% y/u = a/(s+2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Sistema
a     = 3;
num   = a;
den   = [1 2];
ftx   = tf(num,den);

%Dados para simulacao
N     = 500;
dt    = 0.1;
T     = 0:dt:N;
size(T)
u     = ones(1,N/dt+1);

%Simula o sistema para a entrada u degrau unitario
[y,T] = lsim(ftx,u,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% algoritmo do gradiente para identificar parametros
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta(1)     = 0;
erro(1)     = 0;

%Simulação para obter o regressor \phi
ftx1        = tf(1,den);
[phi,T]     = lsim(ftx1,u,T);
gama        = .1;
z = y;

for k = 1:N/dt
    ns2 = phi(k)'*phi(k);
    m2 = 1 + ns2;
    erro(k+1) = (z(k) - theta(k)*phi(k))/m2;
    theta(k+1) = theta(k) + gama*dt*erro(k+1)*phi(k);
end

figure(1)
plot(T,theta)
hold on
plot(T, erro)
legend('Parametro a','Erro')
xlabel('amostras')
ylabel('Parametro identificado')
title('Método do Gradiente')

