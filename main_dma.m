
% MAIN  Inversion of mobility distributions.
% AUTHOR: Timothy Sipkens, 2020-04-11
%=========================================================================%

addpath cmap;

d = logspace(log10(10), log10(1e3), 500)';  % reconstruction points
d_star = logspace(log10(10), log10(1e3), 80)';  % mobility setpoints
d_star2 = exp(log(d_star) - (log(d_star(2)) - log(d_star(2))) ./ 2);  % shifted by 1/2 setpoint

prop_dma = kernel.prop_dma;


A = kernel.gen_dma(d_star, d, prop_dma);


mu_d = 200;
s_d = 1.4;

% Unimodal.
% x0 = normpdf(log10(d), log10(mu_d), log10(s_d));

% Bimodal.
x0 = normpdf(log10(d), log10(mu_d), log10(s_d)) + ...
    0.5 .* normpdf(log10(d), log10(mu_d / 3), log10(s_d / 1.15));

b0 = A * x0;

[b, Lb] = tools.get_noise(b0, 1e2, 1e-6);


figure(1);
plot(d_star, b0, '--', 'Color', [0.85,0.85,0.85]);
hold on;
% plot(d_star, b, 'b.', 'MarkerSize', 10);
stairs(d_star2, b, 'b');  % more representative of TSI display
hold off;
set(gca, 'XScale', 'log');


%%

disp(' ');

%-- Least-squares ---------%
disp('Running least-squares ...');
x_lsq = invert.lsq(Lb*A, Lb*b);
tools.textdone();
disp(' ');


%-- 1st order Tikhonov ----%
disp('Running Tikhonov (1st) ...');
lambda_tk1 = 3.3e1;
[x_tk1, ~, ~, Gpo_inv_tk1] = ...
    invert.tikhonov(Lb*A, Lb*b, lambda_tk1, 1);
Gpo_tk1 = inv(Gpo_inv_tk1);
e.tk1 = (x_tk1 - x0)' * Gpo_inv_tk1 * (x_tk1 - x0);
tools.textdone();
disp(' ');


%-- 2nd order Tikhonov ----%
disp('Running Tikhonov (2nd) ...');
lambda_tk2 = 9e1;
[x_tk2, ~, ~, Gpo_inv_tk2] = ...
    invert.tikhonov(Lb*A, Lb*b, lambda_tk2, 2);
Gpo_tk2 = inv(Gpo_inv_tk2);
e.tk2 = (x_tk2 - x0)' * Gpo_inv_tk2 * (x_tk2 - x0);
tools.textdone();
disp(' ');


%-- Exponential distance --%
disp('Running exponential distance ...');
lambda_ed = 2e0;
ld = 1.2 .* log10(s_d);
[x_ed, ~, ~, Gpo_inv_ed] = ...
    invert.exp_dist(Lb*A, Lb*b, lambda_ed, ld, d);
Gpo_ed = inv(Gpo_inv_ed);
e.ed = (x_ed - x0)' * Gpo_inv_ed * (x_ed - x0);
tools.textdone();
disp(' ');
disp(' ');


e



%%
figure(2);

x_tk = x_tk2;
Gpo_tk = Gpo_tk2;

subplot(1, 2, 1);
tools.plot_ci(d, x_tk, Gpo_tk, x0);
title('Tikhonov');

subplot(1, 2, 2);
tools.plot_ci(d, x_ed, Gpo_ed, x0);
title('Exponential distance');



%%
% Optimize Tikhonov + show Bayes factor.
%{
[a0, a1, a2] = invert.tikhonov_op(Lb*A,Lb*b,[1e-1,1e3],2);

figure(3);
semilogx([a2.lambda], -[a2.B]);
%}

