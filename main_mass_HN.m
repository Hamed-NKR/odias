clc
clear
close all
warning('off')

fdir_aac = 'D:\Hamed\CND\PhD\Publication\Experiment\Aerodynamic_Distribution';
fname_aac = 'aac_dist_feb25';
lbl = {'Lo-Aglom', 'Mod-Colaps', 'Hi-Aglom', 'Ex-Colaps'};
clr = {'#8D493A', '#537188', '#8174A0', '#659287';...
    '#DC8686', '#7EACB5', '#A888B5', '#B1C29E'};
lambda = 0.1;
rd = 2;
n_rsl = 1000;
fdir_dens = 'D:\Hamed\CND\PhD\Publication\Experiment\Effective-Density-Compiled';
fname_dens = 'Effective-Density-Compiled_17-Mar-2025_05-00-44';

load(strcat(fdir_aac, '\', fname_aac, '.mat'), 'dists_aac')
load(strcat(fdir_dens, '\', fname_dens, '.mat'), 'fit4')

n_dat = length(dists_aac);

da = cell(n_dat,1);
dn_dlogda = cell(n_dat,1);
dm = cell(n_dat,1);
rho_eff = cell(n_dat,1);
m = cell(n_dat,1);

f1 = figure(1);
f1.Position = [50, 50, 500, 500];
set(f1, 'color', 'white');

for i = 1 : n_dat
    
    % [dn_dlogda{i}, da{i}] = hn.kde(dists_aac{i}.da,...
    %     dists_aac{i}.dn_dlogda, lambda, rd, n_rsl);
    
    % % Gaussian Process Regression
    % gprMdl = fitrgp(log(dists_aac{i}.da), log(dists_aac{i}.dn_dlogda),...
    %     'KernelFunction', 'squaredexponential');
    % dmin = min(dists_aac{i}.da) / rd;
    % dmax = max(dists_aac{i}.da) * rd;
    % log_da = linspace(log(dmin), log(dmax), n_rsl);
    % da{i} = exp(log_da);
    % dn_dlogda{i} = exp(predict(gprMdl, log_da)); 
    % 
    % plot(da{i}, dn_dlogda{i}, 'Color', hex2rgb(clr{1,i}),...
    %     'LineStyle', '-', 'LineWidth', 1.5);
    hold on
    plot(dists_aac{i}.da, dists_aac{i}.dn_dlogda, 'Color',...
        hex2rgb(clr{2,i}), 'LineStyle', '-.', 'LineWidth', 1.5);

    dm{i} = feval(fit4{i}, dists_aac{i}.da);
    rho_eff{i} = DAT.RHO_EFF(dm{i}, (dists_aac{i}.da)');
    m{i} = (pi/6) * sum(rho_eff{i} .* dm{i}.^3);
    
end

box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.015 0.015], 'XScale', 'log')
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 16)
ylabel('$\mathrm{d}n/\mathrm{dlog}(d_\mathrm{a}) [\#/\mathrm{cm}^3]$',...
    'interpreter', 'latex', 'FontSize', 16)
% xlim([min(cat(2,da{:})), max(cat(2,da{:}))])
% ylim([0,1.05*max(cat(2,dn_dlogda{:}))])


