clear;
close all;
clc;

% import low agglomeration data - repeat 1
load('outputs\26JUL24\26JUL24_rhos.mat')
dm_max_lal_r1 = cat(1, dm_max{:});
rho_max_lal_r1 = cat(1, rho_max{:});

% import high agglomeration data - repeat 1
load('outputs\30JUL24\30JUL24_rhos.mat')
dm_max_hal_r1 = cat(1, dm_max{:});
rho_max_hal_r1 = cat(1, rho_max{:});

% import low agglomeration data - repeat 2
load('outputs\20AUG24_1\20AUG24_1_rhos.mat')
dm_max_lal_r2 = cat(1, dm_max{:});
rho_max_lal_r2 = cat(1, rho_max{:});

% import high agglomeration data - repeat 2
load('outputs\21AUG24\21AUG24_rhos.mat')
dm_max_hal_r2 = cat(1, dm_max{:});
rho_max_hal_r2 = cat(1, rho_max{:});

clear dm_max rho_max % clear redundant variables

% initialize figure
figure();
f1 = gcf;
f1.Position = [100, 100, 500, 500];
set(f1, 'color', 'white');

% plot universal correlation
r0 = (2e4 / 2e0) ^ (1 / (1e4 - 1));
dm0 = 2e0 * ones(1e4, 1) .* r0 .^ (((1:1e4)-1)');
rho0 = 510 * (dm0 / 100) .^ (-0.52);
plt_0 = plot(dm0, rho0, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 2);
hold on

% plot effective density based on mode for low and high agglomeration cases
plt_lal_r1 = scatter(dm_max_lal_r1, rho_max_lal_r1, 15, hex2rgb('#EE4E4E'), 'o');
plt_hal_r1 = scatter(dm_max_hal_r1, rho_max_hal_r1, 15, hex2rgb('#006989'), '^');
plt_lal_r2 = scatter(dm_max_lal_r2, rho_max_lal_r2, 15, hex2rgb('#EE4E4E'), 's');
plt_hal_r2 = scatter(dm_max_hal_r2, rho_max_hal_r2, 15, hex2rgb('#006989'), 'v');

% plot appearance configs
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log',...
    'YminorTick','off')
xlim([50,1000])
xticks(100:100:1000)
xtickangle(90)
ylim([70,700])
yticks(100:100:700)
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 14)
ylabel('$\rho_\mathrm{eff}$ [kg/$\mathrm{m}^{3}$]',...
    'interpreter', 'latex', 'FontSize', 14)
legend(cat(2, plt_0, plt_lal_r1, plt_hal_r1 , plt_lal_r2, plt_hal_r2),...
    cat(2, {'Olfert and Rogak (2019)'},...
    {'Low agglomeration - Repeat 1'}, {'high agglomeration - Repeat 1'},...
    {'Low agglomeration - Repeat 2'}, {'high agglomeration - Repeat 2'}),...
    'interpreter', 'latex', 'FontSize', 14, 'location', 'southwest')
% title(rigstr, 'interpreter', 'latex', 'FontSize', 12)
% subtitle('\textit{1D Tikhonov regularization}', 'interpreter', 'latex',...
%     'FontSize', 12)
exportgraphics(f1, 'outputs\Dens1D_LAL_vs_HAL.png', 'Resolution', 300)

