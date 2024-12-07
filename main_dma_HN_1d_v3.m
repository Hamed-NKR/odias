
clear;
close all;
clc;
warning('off')

addpath cmap; % load cmap library

%% user's inputs to this script %%

tname0 = '13SEP24_SWEEP_ON_FUEL'; % name of excel datasheet for cases to be imported
f_in_add = 'inputs\odias_params_new.xlsx'; % address of folder containing main user info
opts2.correct = 'off'; % flag for whether reiterate the inversion

%% read input files from the user containing details of SMPS data and settings for post-processing %%

opts_load = detectImportOptions(f_in_add);
opts_load = setvartype(opts_load, 'char'); % set type of data to be imported from the excel file
intab = readtable(f_in_add, 'UseExcel', true, 'Sheet', tname0); % read the excel file
fname = intab.fname{1}; % identify name of the file which contains...
% ...information of smps files having raw particle counts
tname = intab.tname{1}; % name of datasheet where the above smps...
% ...information file is located
fname_lbl = intab.fname_lbl{1}; % label for the data to be imported
n_dat = intab.n_dat(1); % number of test cases to be examined
lambda_tk1 = intab.lambda_tk1(1:n_dat); % a parameter that controls smoothing level in Tikhonov method
d_lim_1 = intab.d_lim_1(1:n_dat); % lower limit on mobility diameter (for inversion)
d_lim_2 = intab.d_lim_2(1:n_dat); % upper limit on mobility diameter
clr1 = intab.clr(1:n_dat); % assign color for plot lines
linstl = intab.stl(1:n_dat); % assign line style for plots
dm_rowid = intab.dm_rowid(1:n_dat); % row id in SMPS file that corresponds to mobility setpoints
incld = intab.incld(1:n_dat); % determine whether the class member (data case) undergoes inversion

%% initialize class objects of raw SMPS data %%

dat = DAT(); % initialize the class of SMPS data

% generate a text description for the test
rigstr = DAT.LABELMAKER(fname_lbl, tname);

prop = tfer.prop_dma; % load default SMPS properties

% make placeholders to store mobility setpoints and counts for use in log-scale
dm2 = cell(n_dat, 1);
dN2 = cell(n_dat, 1);

% initialize mobility size distribution figure for all cases
f1 = figure(1);
f1.Position = [100, 100, 1000, 600];
set(f1, 'color', 'white');
t1 = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

plt1 = cell(n_dat,1); % placeholder to store plots for assigning legends

% initialize a plot to check inversion with baseline data at every iteration
f2 = figure(2);
f2.Position = [200, 200, 550, 500];
set(f2, 'color', 'white');

lambda_tk1_0 = lambda_tk1; % store initial inputs for smoothing Tikhonov results
clr2 = {'#FD8A8A', '#F1F7B5', '#A8D1D1', '#9EA1D4'};

% allocate space to TSI inverted size distributions
dist_tsi = struct('fadd', cell(n_dat,1), 'dm', cell(n_dat,1),...
    'dn_dlogdm', cell(n_dat, 1));

% assign class objects to be inverted 
ii = 1 : n_dat;
ii = ii(incld);

tools.textbar([0, n_dat]); % initialize progress textbar


%% iterate through test cases, performing inversion and plotting results %%

for i = ii
    
    obj.f_dir = instab.f_add_raw{pid};
    pid0 = pid;
    while isempty(obj.f_dir) && (pid0 > 0)
        pid0 = pid0 - 1;
        obj.f_dir = instab.f_add_raw{pid0};
    end

    % save the file name, AAC setpoint, AAC and DMA sheath flows...
    % ...and DMA scan indices
    obj.f_name = instab.f_name_raw{pid};

    % import SMPS files and average particle counts
    dat(i) = DAT(i, fname, tname, dm_rowid(i));

    % adjust the description text
    % rigstr{i} = DAT.LABELADJ(rigstr{i});
    
    % reconstruct mobility points
    d = logspace(log10(d_lim_1(i)), log10(d_lim_2(i)), 500)';
    
    % correct sheath flows
    prop.Q_c = dat(i).Qsh_dma;
    prop.Q_m = dat(i).Qsh_dma;

    % build inversion kernel
    A = kernel.gen_smps((dat(i).dm)', d, [], prop);
    
    % assign mean and std of particle counts
    b = (dat(i).dN)';
    Lb = bsxfun(@rdivide, eye(length(b)), (dat(i).sigma)');

    disp(' ');

    % perform inversion using 1st order Tikhonov
    [x_tk1, ~, ~, Gpo_inv_tk1] = ...
        invert.tikhonov(Lb * A, Lb * b, lambda_tk1(i), 1);
    Gpo_tk1 = inv(Gpo_inv_tk1);    

    if strcmp(opts2.correct, 'on') || strcmp(opts2.correct, 'On') ||...
            strcmp(opts2.correct, 'ON')
        
        chk_flag1 = true; % flag for the loop that reiterates the inversion
        chk_flag2 = false;
        
        while chk_flag1
            
            figure(f2) % navigate to the plot that shows interations on inversion

            if chk_flag2
                    [x_tk1, ~, ~, Gpo_inv_tk1] = ...
                        invert.tikhonov(Lb * A, Lb * b, lambda_tk1(i), 1);
                    Gpo_tk1 = inv(Gpo_inv_tk1);
            end
            
            % plot Tikhonov's results based on the last assigned lambda
            plt2_tk = tools.plotci_HN(d, x_tk1, Gpo_tk1, [],...
                hex2rgb(clr2{1}), '-');
            hold on

            % inversion based on Twomey-Markowski
            xi = invert.get_init(Lb * A, Lb * b, d, dat(i).dN);
            x_twomark = invert.twomark(Lb * A, Lb * b, length(xi), xi);
            
            % plot Twomey-Markowski's output (which is closest to TSI's inversion)
            plt2_twomark = tools.plotci_HN(d, x_twomark, [], [],...
                hex2rgb(clr2{4}), '-');
            
            % set the appearance
            set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
                'TickLength', [0.02 0.02], 'XScale', 'log')
            xlim([15,1000])
            xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
                'FontSize', 16)
            legend([plt2_tk, plt2_twomark], {'Tikhonov', 'Twomey-Markowski'},...
                'interpreter', 'latex', 'FontSize', 12, 'Location',...
                'northoutside', 'orientation', 'horizontal')
            
            % prompt to ask user's input on whether inversion is satisfactory
            response2a = questdlg('Does inversion smooting meet the expectations?', ...
	            'Check inversion output', 'Yes', 'No', 'Cancel', 'Yes');

            % handle response
            switch response2a
                case 'Yes'
                    % get out of the loop, plot results and go to the next case
                    chk_flag1 = false;

                case 'No'
                    % generate the prompt message to ask user
                    if ~chk_flag2
                        msg_lambda = strcat('Please enter a refined lambda:');
                        chk_flag2 = true; % flag to redo the inversion
                    else                        
                        msg_lambda = strcat('Please enter a refined lambda',...
                            string(newline), '(initial value:', {' '},...
                            num2str(lambda_tk1_0(i), '%.1e'), ',', {' '},...
                            'recent value:', {' '}, num2str(lambda_tk1(i), '%.1e'),...
                            ')');
                    end

                    % ask user to revise lambda
                    response2b = inputdlg(msg_lambda,...
                        'Refine Tikhonov smoothing', [1 50],...
                        {num2str(lambda_tk1(i))});

                    % terminate the code if user hits cancel
                    if ~isempty(response2b)
                        lambda_tk1(i) = str2double(response2b);
                    else
                        return;
                    end

                case 'Cancel'
                    % terminate the program
                    return
            end

            clf(f2, 'reset')

        end
        
    end

    figure(f1) % navigate to final plot for size distribution
    nexttile(1) % non-log scales    
    
    % display the inverted mobility size distribution based on Tikhonov
    plt1{i} =  tools.plotci_HN(d, x_tk1, Gpo_tk1, [], hex2rgb(clr1{i}),...
        linstl{i});
    hold on
    
    % rigstr{i} = strrep(rigstr{i},'nof','n/f'); % minor label adjustment

    % apply the axis configurations
    if i == n_dat
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
            'TickLength', [0.02 0.02], 'XScale', 'log')    
        xlim([15,1000])
        xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
            'FontSize', 16)
        ylabel('$\mathrm{d}n/\mathrm{dlog}(d_\mathrm{m})$',...
            'interpreter', 'latex', 'FontSize', 16)
    end

    nexttile(2) % log scales

    % remove zero counts for plotting in log scale
    dm2{i} = d(x_tk1 > 0);
    dN2{i} = x_tk1(x_tk1 > 0);
    
    % plot nonzero scales for display in log-log scale
    opts1b.log = 'on';
    tools.plotci_HN(dm2{i}, dN2{i}, Gpo_tk1, [], hex2rgb(clr1{i}),...
        linstl{i}, opts1b)
    hold on

    % apply the log-log scale and other axis configurations
    if i == n_dat
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
            'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')    
        xlim([15,1000])
        ylim([3e3 1.2 * max(cat(1,dN2{:}))])
        xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
            'FontSize', 16)
    end

    tools.textbar([i, n_dat]); % update textbar
end

lgd1 = legend(cat(2, plt1{:}), cat(2, rigstr(:)), 'interpreter', 'latex',...
    'FontSize', 12);
lgd1.Layout.Tile = 'south';
lgd1.NumColumns = 2;


