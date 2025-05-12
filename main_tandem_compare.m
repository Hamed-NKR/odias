clc
clear
close all
warning('off')

%% import data

fdir_in = 'D:\Hamed\CND\PhD\Publication\Experiment\Effective density\MAY25_TANDEM_COMPARE';
fname_in = 'MAY25_TANDEM_COMPARE_08-May-2025_07-17-25';

load(strcat(fdir_in, '\', fname_in, '.mat'), 'dist_odias', 'dist_tsi_inv',...
    'dist_tsi_noninv', 'dat')

%% compare inversion methods

% initialize figure
f1 = figure(1);
f1.Position = [50, 0, 900, 1300];
set(f1, 'color', 'white')
tl1 = tiledlayout(4,3);
tl1.TileSpacing = 'compact';
tl1.Padding = 'compact';

% allocate legend variables
plt1 = cell(6,1);
lgdtxt1 = {'TSI non-corrected', 'TSI inverted',...
     'Twomey-Markowski', 'Exponential distance',...
     '$1^\mathrm{st}$ order Tikhonov', '$2^\mathrm{nd}$ order Tikhonov'};

for i = 1 : 12
    
    nexttile
    hold on

    % TSI without multiple charge correction
    plt{1} = fill([dist_tsi_noninv(i).dm, fliplr(dist_tsi_noninv(i).dm)],...
        [dist_tsi_noninv(i).dn_dlogdm, zeros(size(dist_tsi_noninv(i).dm))],...
        [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.1);

    plt{2} = plot(dist_tsi_inv(i).dm, dist_tsi_inv(i).dn_dlogdm,...
        'LineWidth', 1); % tsi data after inversion
    
    % various tested inversion methods
    plt{3} = plot(dist_odias(i).d, dist_odias(i).x_twomark,...
        'LineWidth', 0.5, 'LineStyle', ':'); % Twomey-Markowski
    plt{4} = plot(dist_odias(i).d, dist_odias(i).x_ed,...
        'LineWidth', 1); % exponential distance
    plt{5} = plot(dist_odias(i).d, dist_odias(i).x_tk1,...
        'LineWidth', 1); % first order Tikhonov
    plt{6} = plot(dist_odias(i).d, dist_odias(i).x_tk2,...
        'LineWidth', 2, 'LineStyle', '-.'); % second order Tikhonov
    
    % set appearance for tiles
    box on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.03 0.03], 'xScale', 'log')
    xlim([70 1000])
    subtitle(strcat('$d_\mathrm{a}$ =', {' '},...
        num2str(dat(i).da, '%.1f'), ' [nm]'),...
        'interpreter', 'latex', 'FontSize', 16);
    
    % print experimental condition
    if i == 2
        title('Lo-Aglom', 'interpreter', 'latex', 'FontSize', 18)
    elseif i == 5
        title('Mod-Colaps', 'interpreter', 'latex', 'FontSize', 18)
    elseif i == 8
        title('Hi-Aglom', 'interpreter', 'latex', 'FontSize', 18)
    elseif i == 11
        title('Ex-Colaps', 'interpreter', 'latex', 'FontSize', 18)
    elseif i == 9
        ylim([0 3e5])
    end

end

% set overall appearance
xlabel(tl1, '$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 24)
ylabel(tl1, '$\mathrm{d}n/\mathrm{dlog}(d_\mathrm{m}) [\#/\mathrm{cm}^3]$',...
    'interpreter', 'latex', 'FontSize', 24)
lgd1 = legend(cat(1, plt1{:}), lgdtxt1, 'interpreter', 'latex',...
    'FontSize', 18, 'Orientation', 'horizontal', 'NumColumns', 3);
lgd1.Layout.Tile = 'north';

%% compare non-inverted data

% initialize figure
f2 = figure(2);
f2.Position = [150, 0, 450, 1300];
set(f2, 'color', 'white')
tl2 = tiledlayout(3,1);
tl2.TileSpacing = 'compact';
tl2.Padding = 'compact';

% allocate legend variables
plt2 = cell(4,1);
lgdtxt2 = {'Lo-Aglom', 'Mod-Colaps', 'Hi-Aglom', 'Ex-Colaps'};

% plot line colors
clr2 = {'#8D493A', '#537188', '#8174A0', '#659287'};

da2 = [100 140 180]; % approximate aac setpoints (for printing)

for i = 1 : 3

    nexttile   
    hold on

    % Lo-Aglom
    plt2{1} = plot(dist_tsi_noninv(i).dm,...
        dist_tsi_noninv(i).dn_dlogdm, 'Color', hex2rgb(clr2{1}),...
        'LineWidth', 2);
    % confidence interval shading (log-safe clipping)
    dm = dist_tsi_noninv(i).dm(:);
    y = dist_tsi_noninv(i).dn_dlogdm(:);
    sigma = dist_tsi_noninv(i).sigma(:);
    ub = y + 2 * sigma;
    lb = max(y - 2 * sigma, 1e2); % avoid log-scale issues
    fill([dm; flipud(dm)], [ub; flipud(lb)],...
        hex2rgb(clr2{1}), 'EdgeColor', 'none', 'FaceAlpha', 0.3);

    % Mod-Colaps
    plt2{2} = plot(dist_tsi_noninv(i+3).dm,...
        dist_tsi_noninv(i+3).dn_dlogdm, 'Color', hex2rgb(clr2{2}),...
        'LineWidth', 2);
    dm = dist_tsi_noninv(i+3).dm(:);
    y = dist_tsi_noninv(i+3).dn_dlogdm(:);
    sigma = dist_tsi_noninv(i+3).sigma(:);
    ub = y + 2 * sigma;
    lb = max(y - 2 * sigma, 1e2);
    fill([dm; flipud(dm)], [ub; flipud(lb)],...
        hex2rgb(clr2{2}), 'EdgeColor', 'none', 'FaceAlpha', 0.3);

    % Hi-Aglom
    plt2{3} = plot(dist_tsi_noninv(i+6).dm,...
        dist_tsi_noninv(i+6).dn_dlogdm, 'Color', hex2rgb(clr2{3}),...
        'LineWidth', 2);
    dm = dist_tsi_noninv(i+6).dm(:);
    y = dist_tsi_noninv(i+6).dn_dlogdm(:);
    sigma = dist_tsi_noninv(i+6).sigma(:);
    ub = y + 2 * sigma;
    lb = max(y - 2 * sigma, 1e2);
    fill([dm; flipud(dm)], [ub; flipud(lb)],...
        hex2rgb(clr2{3}), 'EdgeColor', 'none', 'FaceAlpha', 0.3);

    % Ex-Colaps
    plt2{4} = plot(dist_tsi_noninv(i+9).dm,...
        dist_tsi_noninv(i+9).dn_dlogdm, 'Color', hex2rgb(clr2{4}),...
        'LineWidth', 2);
    dm = dist_tsi_noninv(i+9).dm(:);
    y = dist_tsi_noninv(i+9).dn_dlogdm(:);
    sigma = dist_tsi_noninv(i+9).sigma(:);
    ub = y + 2 * sigma;
    lb = max(y - 2 * sigma, 1e2);
    fill([dm; flipud(dm)], [ub; flipud(lb)],...
        hex2rgb(clr2{4}), 'EdgeColor', 'none', 'FaceAlpha', 0.3);

    % set axes and labels
    box on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16,...
        'TickLength', [0.03 0.03], 'xScale', 'log', 'yScale', 'log')
    xlim([60 1000])
    ylim([1e3 3e6])

    % panel titles
    title(strcat('$d_\mathrm{a}$ =', {' '},...
        num2str(da2(i), '%d'), '$\pm 5$ [nm]'),...
        'interpreter', 'latex', 'FontSize', 18);    
    
end

% set overall appearance
xlabel(tl2, '$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 24)
ylabel(tl2, '$\mathrm{d}n/\mathrm{dlog}(d_\mathrm{m}) [\#/\mathrm{cm}^3]$',...
    'interpreter', 'latex', 'FontSize', 24)
lgd2 = legend(cat(1, plt2{:}), lgdtxt2, 'interpreter', 'latex',...
    'FontSize', 18, 'Orientation', 'horizontal', 'NumColumns', 2);
lgd2.Layout.Tile = 'north';


