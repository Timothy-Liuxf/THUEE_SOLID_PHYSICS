clear; close all; clc;

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

a = 5.43e-10;
m0 = 9.1e-31;
h = 6.63e-34;
N = 9;

% Calculate Fourier of Vn

syms x n;
Vsym(x) = piecewise((-a/4 <= x) & (x <= a/4), 1e-19 .* cos(2*pi/a * x), 0);
Vnsym(n) = 1/a * int(Vsym .* exp(-1j * n .* x * 2*pi/a), x, -a/2, a/2);
Vn = double(vpa(Vnsym(-N : N)));
figure;
hold on;
stem(-N : N, Vn);
title('Fourier');
xlabel('{\itn}');
ylabel('{\itV_n}');
% preprocess

M = floor(N / 2);
selfmap = (-M : M)' - (-M : M) + N + 1;
basemat = Vn(selfmap);
constVal = (h/2/pi)^2 / (2*m0);

nstep = 1000;
k = -pi/a : (2*pi/a / nstep) : pi / a;
lenk = length(k);
eigs = zeros([2 * M + 1, lenk]);
for i = 1 : 1 : lenk
    diagVec = (k(i) - (-M : M) * (2*pi/a)).^2 * constVal;
    A = basemat + diag(diagVec);
    eigs(:, i) = eig(A);
end

figure;
hold on;
kgrid = k / pi * a;
for i = 1 : 1 : 2 * M + 1
    plot(kgrid, eigs(i, :));
end
title('Simple Brillouin picture');
xlabel('k/(pi/{\ita})');
ylabel('{\itE}/J');

midk = round((lenk + 1) / 2);
intval1 = (1 : floor(N / 2) * 2);
for i = 1 : 1 : floor(N / 2)
    intval1(i * 2 - 1) = eigs(i * 2, lenk) - eigs(i * 2 - 1, lenk);
    intval1(i * 2) = eigs(i * 2 + 1, midk) - eigs(i * 2, midk);
end

% Nearly-free electron model

k = -M * pi / a : 2*pi/a/nstep : M * pi / a;
lenk = length(k);
midk = round((lenk + 1) / 2);
kgrid = k / pi * a;

V0 = (Vn((length(Vn) + 1) / 2));
E0 = constVal * k.^2 + V0;
E2 = Vn'.^2 ./ (constVal * (k.^2 - (k + 2*pi/a * (-N:1:N)').^2));
E2(N + 1, :) = 0;
E2 = sum(E2);
E2(isnan(E2)) = 0;
figure;
hold on;
plot(kgrid, E0 + E2);
title('Picture of periodic Brillouin region');
xlabel('k/(pi/{\ita})');
ylabel('{\itE}/J');

E2(abs(E2) > abs(E0)) = 0;
figure;
hold on;
plot(kgrid, E0 + E2);
title('Picture of periodic Brillouin region');
xlabel('k/(pi/{\ita})');
ylabel('{\itE}/J');

idxdelta = pi/2/a/(2*pi/a)*nstep;
E = E0 + E2;
Emid = (length(E) + 1) / 2;
mir = @(idx) Emid - (idx - Emid);
for i = 1 : 1 : M - 1
    ki = i * pi / a;
    [~, idx] = min(abs(k - ki));
    
    [~, maxidx] = max(E(round(idx - idxdelta) : idx - 5));
    maxidx = maxidx + round(idx - idxdelta) - 1;
    E(maxidx:idx) = linspace(E(maxidx), E0(idx) - Vn(i + N + 1), idx - maxidx + 1);

    minidx = idx + (idx - maxidx);
    E(idx:minidx) = linspace(E0(idx) + Vn(i + N + 1), E(minidx), minidx - idx + 1);
    
    [~, idx] = min(abs(k + ki));
    [~, maxidx] = max(E(idx + 5 : round(idx + idxdelta)));
    maxidx = maxidx + idx + 5 - 1;
    E(idx:maxidx) = linspace(E0(idx) - Vn(-i + N + 1), E(maxidx), maxidx - idx + 1);
    
    minidx = idx - (maxidx - idx);
    E(minidx:idx) = linspace(E(minidx), E0(idx) + Vn(-i + N + 1), idx - minidx + 1);
end

figure;
hold on;
plot(kgrid, E);
title('Picture of periodic Brillouin region');
xlabel('k/(pi/{\ita})');
ylabel('{\itE}/J');

tmp = 1 : 8;
intval2(tmp) = 2 * abs(Vn(tmp + N + 1));

figure;
hold on;
stem(intval1);
stem(intval2);
title('Gap');
legend('Simple', 'Period');
xlabel('{\itn}');
ylabel('Gap/J');
