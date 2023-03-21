clear all; clc;
% Parameters - config
a1 = 2;
a2 = -1;
c = 0.5;
%c = 2.0;

theta_1_bar = 2.5; % > abs(a1)
theta_2_bar = 0.5; % > abs(c + a2)


% System matrices
A = [0 1; a1 a2];
b = [0; 1];

x_0 = [-2; 1];
x_i = x_0; %para plot do ponto inicial

x_s = -abs(x_0(1)):abs(x_0(1));
s_plot = -2*x_s;

h = 0.1;

k = 1;

tMax = 100;

nSamples = tMax/h;

while ( k < nSamples)
   s = c*x_0(1) + x_0(2);

   theta_1 = -theta_1_bar * sign(s*x_0(1));
   theta_2 = -theta_2_bar * sign(s*x_0(2));

   u = theta_1*x_0(1) + theta_2*x_0(2);

   x_dot = A*x_0 + b*u;
   x = x_0 + h*x_dot;
   x_0 = x;
   x_1(k) = x(1);
   x_2(k) = x(2);
   k = k + 1;
end

plot(x_1, x_2, 'k-', 'LineWidth', 3);
hold on;
plot(x_s, s_plot, 'r--', 'LineWidth', 3);
hold on;
plot(x_i(1), x_i(2), 'bo', 'LineWidth', 3);
grid on;

