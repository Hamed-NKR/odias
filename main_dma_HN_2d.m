clear;
close all;
clc;
warning('off')

addpath cmap;

%% read user input parameters %%

tname0 = '26JUL24'; % datasheet to be imported
f_in_add = 'inputs\odias_params_2d.xlsx'; 
n_dat0 = 12;

opts = detectImportOptions(f_in_add);
opts = setvartype(opts, 'char');
intab = readtable(f_in_add, 'UseExcel', true, 'Sheet', tname0);
intab = intab(1:n_dat0,:);

fname = intab.fname{1}; % filename to be imported
tname = intab.tname{1}; % tabname that contains data of interest
fname_lbl = intab.fname_lbl{1};
i_dat = intab.i_dat;
lambda_tk1 = intab.lambda_tk1;
d_lim_1 = intab.d_lim_1;
d_lim_2 = intab.d_lim_2;

dat = DAT(); % initialize the class of DMA data

% generate a text description for the test
rigstr = DAT.LABELMAKER(fname_lbl, tname);

% adjust the description text
opts_lbl.bp = 4;
rigstr = DAT.LABELADJ(rigstr{1}, opts_lbl);

% minor label adjustment (for nitrogen premixing part if exists)
if contains(rigstr, 'nof')
    rigstr = strrep(rigstr,'nof','n/f');
end

n_dat = nnz(i_dat);

%% initialize figures %%

% figure for 1d size distributions (in tandem measurements)
figure(1);
f1 = gcf;
f1.Position = [0, 0, ceil(n_dat/16) * 700, 550];
set(f1, 'color', 'white');

%  figure for 1d effective density
figure(2);
f2 = gcf;
f2.Position = [50, 50, ceil(n_dat/16) * 700 + 100, 550];
set(f2, 'color', 'white');

% figure for 2d size distributions
figure(3);
f3 = gcf;
f3.Position = [100, 100, 550, 600];
set(f3, 'color', 'white');

% figure for 2d effective density
figure(4);
f4 = gcf;
f4.Position = [150, 150, 550, 600];
set(f4, 'color', 'white');

% make colormap to distinguish distributions (if input with empty colors)
cm = colormap(turbo); % choose colormap
ii = round(1 + (length(cm) - 1) .* (0.1 : 0.8 / (n_dat - 1) : 0.9));
cm = cm(ii,:); % descretize colors

%% initialize and load variables and parameters for inversion %%

dn_dlogd_max = cell(n_dat,1); % size distribution peaks
dm_max = cell(n_dat,1); % size distribution modes
rho_max = cell(n_dat,1); % effective densities corresponding to modes

dm = cell(n_dat,1); % mobility setpoints
dn_dlogd = cell(n_dat,1); % size distributions
da = cell(n_dat,1); % aerodynamic setpoints

prop = tfer.prop_dma; % load default dma properties

% placeholder to store plots for assigning legends
plt1 = cell(1, n_dat);
plt2 = cell(1, n_dat);
lgd = cell(1, n_dat); % placeholder for legends

j = 0; % index for data to be inverted

tools.textbar([0, n_dat0]); % initialize textbar

%% import and invert tandem measurement data %%

for i = 1 : n_dat0

    % import DMA files and average particle counts
    dat(i) = DAT(i, fname, tname);
        
    if i_dat(i)

        j = j + 1;

        % load colors for size distributions
        if ismatrix(intab.clr)
            if isnan(intab.clr(i)); clr = cm(j,:); end
        else
            clr = hex2rgb(intab.clr{i});
        end

        if isnan(d_lim_1(i)); d_lim_1(i) = min(dat(i).dm); end
        if isnan(d_lim_2(i)); d_lim_2(i) = max(dat(i).dm); end

        % reconstruction points
        dm{i} = logspace(log10(d_lim_1(i)), log10(d_lim_2(i)), 500)';
        
        % correct sheath flows
        prop.Q_c = dat(i).Qsh_dma;
        prop.Q_m = dat(i).Qsh_dma;
    
        % build inversion kernel
        A = kernel.gen_smps((dat(i).dm)', dm{i}, [], prop);
        
        % assign mean and std of particle counts
        b = (dat(i).dN)';
        Lb = bsxfun(@rdivide, eye(length(b)), (dat(i).sigma)');
    
        disp(' ');
    
        % inversion using 1st order Tikhonov
        % disp('Running Tikhonov (1st) ...');
        [dn_dlogd{i}, ~, ~, Gpo_inv_tk1] = ...
            invert.tikhonov(Lb * A, Lb * b, lambda_tk1(i), 1);
        Gpo_tk1 = inv(Gpo_inv_tk1);
        % tools.textdone();
        % disp(' ');
        
        lgd(i) = dat(i).lbl; % record aerodynamic diam. for visualization
        
        % plot mobility size distribution
        figure(f1)
        plt1{i} =  tools.plotci(dm{i}, dn_dlogd{i}, Gpo_tk1, [], clr);
        hold on

        % find modes of distribution
        [dn_dlogd_max{i}, dm_max{i}] = findpeaks(dn_dlogd{i}, dm{i});

        % remove noises from the modes
        dm_max{i}(dn_dlogd_max{i} / max(dn_dlogd_max{i}) < 0.1) = [];
        dn_dlogd_max{i}(dn_dlogd_max{i} / max(dn_dlogd_max{i}) < 0.1) = [];
        
        % remove out of range peaks
        dn_dlogd_max{i}(dm_max{i} > 1000) = [];
        dm_max{i}(dm_max{i} > 1000) = [];

        % find 1d effective densities based on modes
        rho_max{i} = DAT.RHO_EFF(dm_max{i}, dat(i).da);
        
        % reformat aerodynamic setpoint
        da{i} = repmat(dat(i).da, length(dm{i}), 1);
        
        % plot 1d effective density
        figure(f2)
        plt2{i} = scatter(dm_max{i}, rho_max{i}, 15, clr, 'o');
        hold on

    end
    
    tools.textbar([i, n_dat0]); % update textbar
    
end


%% plot appearance configs for 1d size distributions %%

figure(f1);
ax1 = gca;

set(ax1, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log')

xlim([30,1000])

xlabel(ax1, '$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 14)
ylabel(ax1, '$\mathrm{d}n/\mathrm{dlog}(d_\mathrm{m})$',...
    'interpreter', 'latex', 'FontSize', 14)

legend(cat(2, plt1{i_dat}), cat(2, lgd(i_dat)), 'interpreter', 'latex',...
    'FontSize', 14, 'location', 'eastoutside', 'NumColumns', ceil(n_dat/16))

title(ax1, rigstr, 'interpreter', 'latex', 'FontSize', 12)
subtitle(ax1, '\textit{1D Tikhonov regularization}', 'interpreter', 'latex',...
    'FontSize', 12)

exportgraphics(f1, strcat('outputs\', tname0, '_TandemDist.png'),...
    'Resolution', 300)

%% configs for 1d effective density plot %%

figure(f2);

% plot universal correlation
r0 = (2e4 / 2e0) ^ (1 / (1e4 - 1));
dm0 = 2e0 * ones(1e4, 1) .* r0 .^ (((1:1e4)-1)');
rho0 = 510 * (dm0 / 100) .^ (-0.52);
plt2_0 = plot(dm0, rho0, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 2);

ax2 = gca;
box on

set(ax2, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log',...
    'YminorTick','off')

xlim([50,1000])
xticks(100:100:1000)
xtickangle(90)
ylim([70,700])
yticks(100:100:700)

xlabel(ax2, '$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 14)
ylabel(ax2, '$\rho_\mathrm{eff}$ [kg/$\mathrm{m}^{3}$]',...
    'interpreter', 'latex', 'FontSize', 14)

legend(cat(2, plt2{i_dat}, plt2_0), cat(2, lgd(i_dat), {'Olfert and Rogak (2019)'}),...
    'interpreter', 'latex', 'FontSize', 14, 'location', 'eastoutside',...
    'NumColumns', ceil(n_dat/16))

title(ax2, rigstr, 'interpreter', 'latex', 'FontSize', 12)
subtitle(ax2, '\textit{1D Tikhonov regularization}', 'interpreter', 'latex',...
    'FontSize', 12)

exportgraphics(f2, strcat('outputs\', tname0, '_Dens1D.png'),...
    'Resolution', 300)

%% generate 2d size distribution heatmap %%

figure(f3);

dm2 = cat(1, dm{:}); % 2d map for mobility setpoints
dn2_dlogd = cat(1, dn_dlogd{:}); % 2d map for size distributions
da2 = cat(1, da{:}); % aerodynamic diameters in 2d space

cm2 = colormap(hot); % choose colormap
% cm2 = flip(cm2,1); % flip the order of colors if needed

% assign size distributions to colormap
dn2_dlogd_hat = round(rescale(log(1 + dn2_dlogd), 1, 256));
dn2_dlogd_hat(isnan(dn2_dlogd_hat)) = 1;

% plot a scattered heatmap of mobility vs. aerodynamic diameter...
    % ...colored based on size distribution
scatter(dm2, da2, 50, cm2(dn2_dlogd_hat,:), 's', 'filled');

ax3 = gca;
box on

set(ax3, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')

xlim([30, 1000])
ylim([30, 400])

% print colorbar
cb3 = colorbar;
cb3.Label.String = '$\mathrm{log}\left(1 + ({\mathrm{d}^2}n/{\mathrm{dlog}(d_\mathrm{m})}{\mathrm{dlog}(d_\mathrm{a})})/({\mathrm{d}^2}n_\mathrm{max}//{\mathrm{dlog}(d_\mathrm{m})}{\mathrm{dlog}(d_\mathrm{a})})\right)$';
cb3.Label.Interpreter  = 'latex';
% cb3.Label.Rotation = 360;
cb3.TickLabelInterpreter  = 'latex';
cb3.FontSize = 12;
cb3.Label.FontSize = 12;
% cb3.Limits = [1.0 1.6];
% cb3.Ticks = 1.0:0.1:1.6;
cb3.TickLength = 0.02;
cb3.Location = 'southoutside';
% cb3pos = get(cb3, 'Position');
% cb3.Label.Position = [cb3pos(1) - 0.75 , cb3pos(2) + 1.705];

xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex', 'FontSize', 16)
ylabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 16)

title(ax3, rigstr, 'interpreter', 'latex', 'FontSize', 12)
subtitle(ax3, '\textit{1D Tikhonov regularization}', 'interpreter', 'latex',...
    'FontSize', 12)

exportgraphics(f3, strcat('outputs\', tname0, '_Dist2d.png'),...
    'Resolution', 300)

%% generate 2d effective density heatmap %%

% calculate effective density distribution in 2d space
rho2 = DAT.RHO_EFF(dm2, da2);

figure(f4);
ax4 = gca;

% scattered heatmap for effective density vs. mobility size
scatter(dm2, rho2, 50, cm2(dn2_dlogd_hat,:), 's', 'filled');

hold on

% plot universal correlation
plt4_0 = plot(dm0, rho0, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 2);

box on
set(ax4, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')

xlim([30, 1000])
ylim([20, 2000])

colormap(hot);

% print colorbar
cb4 = colorbar;
cb4.Label.String = '$\mathrm{log}\left(1 + ({\mathrm{d}^2}n/{\mathrm{dlog}(d_\mathrm{m})}{\mathrm{dlog}(d_\mathrm{a})})/({\mathrm{d}^2}n_\mathrm{max}//{\mathrm{dlog}(d_\mathrm{m})}{\mathrm{dlog}(d_\mathrm{a})})\right)$';
cb4.Label.Interpreter  = 'latex';
% cb.Label.Rotation = 360;
cb4.TickLabelInterpreter  = 'latex';
cb4.FontSize = 12;
cb4.Label.FontSize = 12;
% cb4.Limits = [1.0 1.6];
% cb4.Ticks = 1.0:0.1:1.6;
cb4.TickLength = 0.02;
cb4.Location = 'southoutside';
% cb4pos = get(cb4, 'Position');
% cb4.Label.Position = [cb4pos(1) - 0.75 , cb4pos(2) + 1.705];

xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex', 'FontSize', 16)
ylabel('$\rho_\mathrm{eff}$ [kg/$\mathrm{m}^{3}$]', 'interpreter', 'latex', 'FontSize', 16)

legend(plt4_0, 'Olfert and Rogak (2019)', 'interpreter', 'latex',...
    'FontSize', 12, 'location', 'southwest')

title(ax4, rigstr, 'interpreter', 'latex', 'FontSize', 12)
subtitle(ax4, '\textit{1D Tikhonov regularization}', 'interpreter', 'latex',...
    'FontSize', 12)

exportgraphics(f4, strcat('outputs\', tname0, '_Dens2d.png'),...
    'Resolution', 300)

%% copy 1d size distributions to show in log-log scale %%

figure(f1);
lgd1 = legend('show');
f1_props = get(f1);
f5 = figure('Position', f1_props.Position, ...
    'Color', f1_props.Color, ...
    'Colormap', f1_props.Colormap, ...
    'Renderer', f1_props.Renderer);
ax5 = copyobj([lgd1, ax1], f5);
ttl1 = get(ax1, 'Title');
ttl5 = get(ax5, 'Title');
ttl5 = ttl5{2};
set(ttl5, 'FontSize', get(ttl1, 'FontSize'));
sub1 = get(ax1, 'Subtitle');
sub5 = get(ax5(2), 'Subtitle');
set(sub5, 'FontSize', get(sub1, 'FontSize'));
set(gca, 'YScale', 'log')
ylim([1, inf])

exportgraphics(f5, strcat('outputs\', tname0, '_TandemDist_Log.png'),...
    'Resolution', 300)

%% export workspace data %%

save(strcat('outputs\', tname0)) % save all data

% save only 1d effective density data
save(strcat('outputs\', tname0, '_rhos'), 'rho_max', 'dm_max')
