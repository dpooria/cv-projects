clc;
clear;
%-- import data
dat = dlmread('fitinput.txt');
x_ob = dat(:, 1); y_ob = dat(:, 2);
%-- find a_best and H_best with matlab fit function
y_t = @(a,H,x)(a * x.^H);
ft = fit(x_ob, y_ob, y_t);
a_best = ft.a; H_best = ft.H;
%-----------------------
syms a H x y;
vay = var(y_ob);
X2 = (y - a * x^H)^2;
X2_min = (y - a_best * x^H_best)^2; % constant
DX2 = (X2 - X2_min)/vay;
% F_ij = 0.5 * d^2(DX^2)/dthi dthj |_th = th_best
F11 = 0.5 * subs(diff(DX2, a, 2), [a, H], [a_best, H_best]);
F12 = 0.5 * subs(diff(diff(DX2, a), H), [a,H], [a_best, H_best]);
F21 = 0.5 * subs(diff(diff(DX2, H), a), [a, H], [a_best, H_best]);
F22 = 0.5 * subs(diff(DX2, H, 2), [a, H], [a_best, H_best]);
F = zeros(2);
N = length(x_ob);
for i = 1:N
    f11 = subs(F11, [x, y], [x_ob(i), y_ob(i)]);
    f12 = subs(F12, [x, y], [x_ob(i), y_ob(i)]);
    f21 = subs(F21, [x, y], [x_ob(i), y_ob(i)]);
    f22 = subs(F22, [x, y], [x_ob(i), y_ob(i)]);
    % F ~ < [F11, F12; F21, F22] >
    F = F + double([f11, f12; f21, f22]);
end
F = F / N;
disp('Fisher Matrix:');
disp(F);
%-- inverse of Fisher Matrix gives the systematic errors of a and H
Ct = abs(inv(F));
disp('Error_a:');
disp(Ct(1,1));
disp('Error_H:');
disp(Ct(2,2));
disp('Error_aH:');
disp(Ct(1,2));
disp('Error_Ha:');
disp(Ct(2,1));
%-- calculate DX^2 with elliptical equation
% p -> pearson coefficient
p = Ct(1,2) / sqrt(Ct(1,1) * Ct(2,2));
n = 1000;
da = linspace(-10, 10, n) - a_best;
dH = linspace(-10, 10, n) - H_best;
DX2_el = zeros(n);
for i = 1:n
    for j = 1:n
        DX2_el(i,j) = (da(i)^2)/(Ct(1,1)*(1-p^2)) + ...
            (dH(j)^2)/(Ct(2,2)*(1-p^2)) - ...
            (2*p/(1-p^2)) * da(i) * dH(i) / sqrt(Ct(1,1) * Ct(2,2));
    end
end
% calculate DX^2 with Fisher formalism
dth = [da; dH];
DX2_f = dth' * F * dth;
% plot DX^2
fig1 =  figure();
contour(da, dH, DX2_el);
grid on;
xlabel('\deltaa');
ylabel('\deltaH');
title('\DeltaX^2');
colorbar();
saveas(fig1, 'DX2_el.png');
fig2 = figure();
contour(da, dH, DX2_f);
grid on;
xlabel('\deltaa');
ylabel('\deltaH');
title('\DeltaX^2');
colorbar();
saveas(fig2, 'DX2_f.png');
