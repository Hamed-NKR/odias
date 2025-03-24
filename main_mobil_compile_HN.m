clc
clear
close all
warning('off')

fdir = {'D:\Hamed\CND\PhD\Publication\Experiment\Mobility_Distribution\ET017_NIT013_AIR16_DIST',...
    'D:\Hamed\CND\PhD\Publication\Experiment\Mobility_Distribution\MOBIL_DIST_FEB25'};
fname = {'ET017_NIT013_AIR16_DIST_27-Dec-2024_23-29-11',...
    'MOBIL_DIST_FEB25_15-Mar-2025_21-31-39'};
ii = {1:7,8:16,[],[];1:5,6,7:10,11:13};
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
dists_dma = cell(n_dat,1);

f1 = figure(1);
f1.Position = [100, 100, 500, 500];
set(f1, 'color', 'white');
plt1 = cell(n_dat,1);

for i = 1 : n_dat
    
    dists_dma{i}.dm0 = cat(2,dist_odias_1(ii{1,i}).d, dist_odias_2(ii{2,i}).d);
    dists_dma{i}.dn0_dlogdm = cat(2,dist_odias_1(ii{1,i}).x, dist_odias_2(ii{2,i}).x);
    
    n_dm0 = size(dists_dma{i}.dm0);
    [dists_dma{i}.dn_dlogdm, dists_dma{i}.dm, dists_dma{i}.ci_dn_dlogdm] = ...
        hn.avg_dist_2(mat2cell(dists_dma{i}.dm0, n_dm0(1), ones(n_dm0(2),1)),...
        mat2cell(dists_dma{i}.dn0_dlogdm, n_dm0(1), ones(n_dm0(2),1)));

    % compute the bounds
    dists_dma{i}.ub_dn_dlogdm = dists_dma{i}.dn_dlogdm +...
        dists_dma{i}.ci_dn_dlogdm;
    dists_dma{i}.lb_dn_dlogdm = dists_dma{i}.dn_dlogdm -...
        dists_dma{i}.ci_dn_dlogdm;
    
    iii = dists_dma{i}.dn_dlogdm < 0;

    if nnz(iii)
        dists_dma{i}.dm(iii) = [];
        dists_dma{i}.dn_dlogdm(iii) = [];
        dists_dma{i}.ci_dn_dlogdm(iii) = [];
        dists_dma{i}.ub_dn_dlogdm(iii) = [];
        dists_dma{i}.lb_dn_dlogdm(iii) = [];
    end

    plt1{i} = plot(dists_dma{i}.dm, dists_dma{i}.dn_dlogdm,...
        'Color', hex2rgb(clr{1,i}), 'LineWidth', 1.5);
    hold on
    fill([dists_dma{i}.dm, fliplr(dists_dma{i}.dm)],...
        [dists_dma{i}.ub_dn_dlogdm, fliplr(dists_dma{i}.lb_dn_dlogdm)],...
        hex2rgb(clr{2,i}), 'EdgeColor', 'none', 'FaceAlpha', 0.5);

    % total concentration
    dists_dma{i}.n_tot = trapz(log10(dists_dma{i}.dm),...
        dists_dma{i}.dn_dlogdm);
    
    % weights for calculating geometric mean
    w = dists_dma{i}.dn_dlogdm / sum(dists_dma{i}.dn_dlogdm);
    
    % geometric mean
    dists_dma{i}.gm_dm = 10^(sum(w .* log10(dists_dma{i}.dm)));
    
    % geometric standard deviation
    dists_dma{i}.gsd_dm = 10^(sqrt(sum(w .* (log10(dists_dma{i}.dm) -...
        log10(dists_dma{i}.gm_dm)).^2)));    
    
end

% set apprearances
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.015 0.015], 'XScale', 'log')
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 16)
ylabel('$\mathrm{d}n/\mathrm{dlog}(d_\mathrm{m}) [\#/\mathrm{cm}^3]$',...
    'interpreter', 'latex', 'FontSize', 16)
xlim([25,1000])
ylim([0,1.4e8])
legend(cat(2, plt1{:}), cat(2, lbl(:)), 'interpreter', 'latex',...
    'FontSize', 12, 'Location', 'northeast');

