
clear;
close all;
clc;
warning('off')

addpath cmap;

%%% load data %%%
 
dat = DAT(); % initialize the class of DMA data

fname = '20JUL24\DMA_PARAMS'; % filename to be imported
tname = '20JUL24'; % tabname that contains data of interest
fname_props = '20JUL24\PFA_AB_PROPS';
n_dat = 7;

% generate a text description for the test
rigstr = DAT.LABELMAKER(fname_props, tname);

% import DMA files, average and assign mobility data to AAC setpoints
for i = 1 : n_dat
    dat(i) = DAT(i, fname, tname);

    % adjust the description text
    rigstr{i} = DAT.LABELADJ(rigstr{i});
end

%%%%%%

prop = tfer.prop_dma;
lambda_tk1 = [5*1e-3, 5*1e-3, 5*1e-3, 5*1e-3, 5*1e-3, 5*1e-3, 5*1e-3];
clr = {'#1b85b8' , '#5a5255', '#559e83', '#ae5a41', '#c3cb71',...
    '#987D9A', '#F19ED2', '#C7C8CC'};
d_lim = [15, 1.5e3; 15, 1.5e3; 15, 1.5e3; 15, 1.5e3; 15, 1.5e3;...
    15, 1.5e3; 15, 1.5e3];

figure;
h = gcf;
h.Position = [0, 0, 900, 500];
set(h, 'color', 'white');

plt = cell(n_dat,1);

for i = 1 : n_dat
    
    d = logspace(log10(d_lim_1(i)), log10(d_lim_2(i)), 500)';  % reconstruction points

    prop.Q_c = dat(i).Qsh_dma;
    prop.Q_m = dat(i).Qsh_dma;
    
    A = kernel.gen_smps((dat(i).dm)', d, [], prop);

    b = (dat(i).dN)';
    Lb = bsxfun(@rdivide, eye(length(b)), (dat(i).sigma)');

    disp(' ');

    %-- 1st order Tikhonov ----%
    disp('Running Tikhonov (1st) ...');
    
    [x_tk1, ~, ~, Gpo_inv_tk1] = ...
        invert.tikhonov(Lb * A, Lb * b, lambda_tk1(i), 1);
    Gpo_tk1 = inv(Gpo_inv_tk1);
    tools.textdone();
    disp(' ');

    plt{i} =  tools.plotci(d, x_tk1, Gpo_tk1, [], hex2rgb(clr{i}));
    hold on
    
    rigstr{i} = strrep(rigstr{i},'nof','n/f');
end

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15,...
    'TickLength', [0.02 0.02])
xlim([15,1500])
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 18)
ylabel('$\mathrm{d}N/Q_\mathrm{DMA}$ [$\mathrm{\#/cm\hat{A}^3}$]',...
    'interpreter', 'latex', 'FontSize', 18)
legend(cat(2, plt{:}), cat(2, rigstr{:}), 'interpreter', 'latex',...
    'Location', 'eastoutside', 'NumColumns', 1, 'FontSize', 15)


