clc
clear
close all
warning('off')

fdir = {'D:\Hamed\CND\PhD\Publication\Experiment\Mobility_Distribution\ET017_NIT013_AIR16_DIST',...
    'D:\Hamed\CND\PhD\Publication\Experiment\Mobility_Distribution\MOBIL_DIST_FEB25'};
fname = {'ET017_NIT013_AIR16_DIST_27-Dec-2024_23-29-11',...
    'MOBIL_DIST_FEB25_15-Mar-2025_21-31-39'};
ii = {1:7,8:16,[],[];1:5,6,7:10,11:13};
cm = {'spring','summer','autumn','winter'};
lbl = {'Lo-Aglom', 'Mod-Colaps', 'Hi-Aglom', 'Ex-Colaps'};
clr = {'#8D493A', '#537188', '#8174A0', '#659287';...
    '#DC8686', '#7EACB5', '#A888B5', '#B1C29E'};

for j = 1:length(fname)
    
    load(strcat(fdir{j}, '\', fname{j}, '.mat'), 'dist_odias')
    
    varname = who;
    jj = find(strcmp(varname, 'dist_odias'));
    newVarName = [varname{jj} '_' num2str(j)];
    eval([newVarName ' = ' varname{jj} ';']);
    eval(['clear ', varname{jj}])

end

clear varname newVarName

n_dat = length(lbl);
dm0 = cell(n_dat,1);
dn0_dlogdm = cell(n_dat,1);
dm = cell(n_dat,1);
dn_dlogdm = cell(n_dat,1);
ub_dn_dlogdm = cell(n_dat,1);
lb_dn_dlogdm = cell(n_dat,1);

f1 = figure(1);
f1.Position = [100, 100, 600, 600];
set(f1, 'color', 'white');
plt1 = cell(n_dat,1);

for i = 1 : n_dat
    
    dm0{i} = cat(2,dist_odias_1(ii{1,i}).d, dist_odias_2(ii{2,i}).d);
    dn0_dlogdm{i} = cat(2,dist_odias_1(ii{1,i}).x, dist_odias_2(ii{2,i}).x);
    
    n_dm0 = size(dm0{i});
    [dn_dlogdm{i}, dm{i}, ci_dn_dlogdm{i}] = ...
        hn.avg_dist_2(mat2cell(dm0{i}, n_dm0(1), ones(n_dm0(2),1)),...
        mat2cell(dn0_dlogdm{i}, n_dm0(1), ones(n_dm0(2),1)));

    % compute the bounds
    ub_dn_dlogdm{i} = dn_dlogdm{i} + ci_dn_dlogdm{i};
    lb_dn_dlogdm{i} = dn_dlogdm{i} - ci_dn_dlogdm{i};
    
    % iii = (dn_dlogdm{i} <= 0) | (ub_dn_dlogdm{i} <= 0) |...
    %     (lb_dn_dlogdm{i} <= 0);
    % dm{i}(iii) = [];
    % dn_dlogdm{i}(iii) = [];
    % ub_dn_dlogdm{i}(iii) = [];
    % lb_dn_dlogdm{i}(iii) = [];

    plt1{i} = plot(dm{i}, dn_dlogdm{i}, 'Color', hex2rgb(clr{1,i}),...
        'LineWidth', 1.5);
    hold on
    fill([dm{i}, fliplr(dm{i})],...
        [ub_dn_dlogdm{i}, fliplr(lb_dn_dlogdm{i})],...
        hex2rgb(clr{2,i}), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    
end

% set apprearances
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.015 0.015], 'XScale', 'log')
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 16)
ylabel('$\mathrm{d}n/\mathrm{dlog}(d_\mathrm{m}) [\#/\mathrm{cm}^3]$',...
    'interpreter', 'latex', 'FontSize', 16)
xlim([25,1000])
ylim([0,1.025*max(cat(2,ub_dn_dlogdm{:}))])
legend(cat(2, plt1{:}), cat(2, lbl(:)), 'interpreter', 'latex',...
    'FontSize', 12, 'Location', 'northoutside', 'Orientation', 'horizontal');

