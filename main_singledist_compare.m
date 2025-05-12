clc
clear
close all
warning('off')

%% import data and initialize

% import AAC data
fdir_aac = 'D:\Hamed\CND\PhD\Publication\Experiment\Aerodynamic_Distribution';
fname_aac = 'aac_dist_feb25';
load(strcat(fdir_aac, '\', fname_aac, '.mat'), 'dists_aac')

% import SMPS data
fdir_smps = 'D:\Hamed\CND\PhD\Publication\Experiment\Mobility_Distribution\MOBIL_DIST_COMPILED_MAY25';
fname_smps = 'MOBIL_DIST_FEB25_12-May-2025_13-24-39';
load(strcat(fdir_smps, '\', fname_smps, '.mat'), 'dists_dma')

% initialize figure
f1 = figure(1);
f1.Position = [50, 50, 1000, 550];
set(f1, 'color', 'white')
tl1 = tiledlayout(1,2);
tl1.TileSpacing = 'compact';
tl1.Padding = 'compact';

% allocate legend variables
plt1 = cell(4,1);
lgdtxt1 = {'Lo-Aglom', 'Mod-Colaps', 'Hi-Aglom', 'Ex-Colaps'};
clr1 = {'#8D493A', '#537188', '#8174A0', '#659287'};

%% compare one-dimensional "AAC" distributions

nexttile(1)

for i = 1 : 4
    
    % average distribution
    plt1{i} = plot(dists_aac{i}.da, dists_aac{i}.dn_dlogda, 'Color',...
        hex2rgb(clr1{i}), 'LineWidth', 1.5);
    hold on

    % confidence intervals
    ub_dn_dlogda = dists_aac{i}.dn_dlogda + dists_aac{i}.ci_dn_dlogd;
    lb_dn_dlogda = max(dists_aac{i}.dn_dlogda - dists_aac{i}.ci_dn_dlogd,...
        1e2);
    fill([dists_aac{i}.da, fliplr(dists_aac{i}.da)],...
        [ub_dn_dlogda, fliplr(lb_dn_dlogda)],...
        hex2rgb(clr1{i}), 'EdgeColor', 'none', 'FaceAlpha', 0.4);    
end

% set figure appearances
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16,...
    'TickLength', [0.02 0.02], 'xScale', 'log', 'yScale', 'log')
xlim([35,800])
xticks([50 100 200 400 800])
xtickangle(90)
xlabel('$d_\mathrm{ae}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 24)
ylim([1e4, 1.5e8])
ylabel('$\mathrm{d}n/\mathrm{dlog}(d_\mathrm{ae}) [\#/\mathrm{cm}^3]$',...
    'interpreter', 'latex', 'FontSize', 24)
lgd1 = legend(cat(1, plt1{:}), lgdtxt1, 'interpreter', 'latex',...
    'FontSize', 20, 'Orientation', 'horizontal');
lgd1.Layout.Tile = 'south';

%% compare one-dimensional "SMPS" distributions

nexttile(2)

for i = 1 : 4
    
    % average distribution
    plot(dists_dma{i}.dm, dists_dma{i}.dn_dlogdm, 'Color',...
        hex2rgb(clr1{i}), 'LineWidth', 1.5)
    hold on

    % confidence intervals
    fill([dists_dma{i}.dm, fliplr(dists_dma{i}.dm)],...
        [dists_dma{i}.ub_dn_dlogdm,...
        fliplr(max(dists_dma{i}.lb_dn_dlogdm, 1e2))],...
        hex2rgb(clr1{i}), 'EdgeColor', 'none', 'FaceAlpha', 0.4); 

end

% set figure appearances
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16,...
    'TickLength', [0.02 0.02], 'xScale', 'log', 'yScale', 'log')
xlim([18,1000])
xticks([20 40 60 80 100 200 400 600 800 1000])
xtickangle(90)
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 24)
ylim([1e3, 1.5e8])
ylabel('$\mathrm{d}n/\mathrm{dlog}(d_\mathrm{m}) [\#/\mathrm{cm}^3]$',...
    'interpreter', 'latex', 'FontSize', 24)

