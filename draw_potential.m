clear; close; clc;

% constants

a = 5.43e-10;
m0 = 9.1e-31;
h = 6.63e-34;

Vfunc = @(x) (2*pi/a * x >= -pi/2 & 2*pi/a * x <= pi/2) .* 1e-19 .* cos(2*pi/a * x);
figure;
hold on;
fplot(@(x) Vfunc(a/(2*pi)*x) * 1e19, [-pi, pi]);
title('Potential Energy Distribution');
xlabel('2\pi{\itx}/{\ita}');
ylabel('{\itV}/(10^{-19}J)')
