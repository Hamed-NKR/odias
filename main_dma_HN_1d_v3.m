
clear;
close all;
clc;
warning('off')

addpath cmap; % load cmap library

%% user's inputs to this script %%

% assign name of excel datasheet for cases to be imported
% tname0 = '13SEP24_SWEEP_ON_FUEL'; % parametric study on effect of...
    % % ...fuel flow rate (i.e. change of equivalence ratio)
% tname0 = '20JUL24_SWEEP_ON_NITROGEN'; % effect of change in premixed...
    % % ...Nitrogen to fuel mass ratio 
% tname0 = '23JUL24_SWEEP_ON_AIR'; % effect of change in air flow rate...
    % % ...while equivalence ratio is constant
% tname0 = 'ET017_NIT013_AIR16_DIST'; % size distributions for...
%     % ...favorite burner condition
tname0 = '19AUG24_LOWAGGLOM_TANDEM'; % tandem measurements for low...
    % ...agglomeration case on 19 Aug. 2024

% address of folder containing main user info
f_in_add = 'inputs\odias_params_new.xlsx';

opts2.correct = 'on'; % flag for whether reiterate the inversion

%% read input files from the user containing details of SMPS data and settings for post-processing %%

opts_load = detectImportOptions(f_in_add);
opts_load = setvartype(opts_load, 'char'); % set type of data to be...
    % ...imported from the user-info excel file
intab = readtable(f_in_add, 'UseExcel', true, 'Sheet', tname0); % read the excel file

fname = intab.fname{1}; % identify name of the file which contains...
    % ...information of smps files having raw particle counts
tname = intab.tname{1}; % name of datasheet where the above smps...
    % ...information file is located
fname_lbl = intab.fname_lbl{1}; % label for the data to be imported
n_dat = intab.n_dat(1); % number of test cases to be examined

ii = 1 : n_dat; % initially all data are read

lambda_tk1 = intab.lambda_tk1(ii); % a parameter that controls...
    % ...smoothing level in first-order Tikhonov method
lambda_tk2 = intab.lambda_tk2(ii); % smoothing in 2nd-order Tikhonov
lambda_ed = intab.lambda_ed(ii); % smoothing by exponential distance
d_lim_1 = intab.d_lim_1(ii); % lower limit on mobility diameter (for inversion)
d_lim_2 = intab.d_lim_2(ii); % upper limit on mobility diameter
clr1 = intab.clr(ii); % assign color for plot lines
linstl = intab.stl(ii); % assign line style for plots
dm_rowid = intab.dm_rowid(ii); % row id in SMPS file that...
    % ...corresponds to mobility setpoints
incld = intab.incld(ii); % determine whether the class member...
    % ...(data case) undergoes inversion

% assign class objects to be inverted 
ii = ii(incld);

% group similar data (if input by user) for possible later averaging and...
    % ...colorcoding
if ismember('group', fieldnames(intab))

    grp = intab.group(ii); % load grouping indices

    % find unique group indices
    uniq_grp = unique(grp(ii));

    n_grp = numel(uniq_grp); % number of unique groups
    
    % initialize placeholders for indices of datasets corresponding to...
        % ...each group
    ind_grp = cell(n_grp, 1);
    
    % find indices of members for each group
    for j = 1 : n_grp
        ind_grp{j} = find(grp(ii) == uniq_grp(j));
    end

end

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

% store initial smoothing and size parameters for Tikhonov inversion
lambda_tk1_0 = lambda_tk1;
lambda_tk2_0 = lambda_tk2;
lambda_ed_0 = lambda_ed;
d_lim_1_0 = d_lim_1;
d_lim_2_0 = d_lim_2;

% colors for plot that checks inversion smoothing
clr2 = {'#295F98', '#659287', '#A888B5', '#C96868', '#EF9C66' ,'#FCDC94',...
    '#C8CFA0'};

% colors for average size distribution plots
clr3 = {{'#8D493A', '#DC8686'}, {'#537188', '#7EACB5'}};

% allocate space to TSI inverted size distributions
dist_tsi_inv = struct('fadd', cell(n_dat,1), 'tab0', cell(n_dat,1),...
    'dm', cell(n_dat,1), 'dn0_dlogdm', cell(n_dat, 1),...
    'dn_dlogdm', cell(n_dat, 1), 'sigma', cell(n_dat, 1),...
    'n0_tot', cell(n_dat, 1), 'n_tot', zeros(n_dat,1),...
    'A_tot', zeros(n_dat,1));
dist_tsi_noninv = struct('fadd', cell(n_dat,1), 'tab0', cell(n_dat,1),...
    'dm', cell(n_dat,1), 'dn0_dlogdm', cell(n_dat, 1),...
    'dn_dlogdm', cell(n_dat, 1), 'sigma', cell(n_dat, 1),...
    'n0_tot', cell(n_dat, 1), 'n_tot', zeros(n_dat,1),...
    'A_tot', zeros(n_dat,1));

% allocate space to Tikhonov inversion results
dist_odias = struct('d', cell(n_dat,1), 'x', cell(n_dat,1),...
    'Gpo', cell(n_dat, 1), 'x_twomark', cell(n_dat,1),...
    'x_tk1', cell(n_dat,1), 'Gpo_tk1', cell(n_dat, 1),...
    'x_tk2', cell(n_dat,1), 'Gpo_tk2', cell(n_dat, 1),...
    'x_ed', cell(n_dat,1), 'Gpo_ed', cell(n_dat, 1),...
    'n_tot', zeros(n_dat,1), 'd_mode', cell(n_dat, 1),...
    'd_gm', zeros(n_dat, 1), 'sigma_g', cell(n_dat, 1),...
    'skw', zeros(n_dat, 1), 'krts', zeros(n_dat, 1));

% define a multiplier to correct for conversion issues of inversion algorithm
rN = zeros(n_dat,4); % first column: Tikhonov 1st order,...
    % ...second column: Tikhonov 2nd order,...
    % ...third column: Tikhonov 2nd order and two-setp,...
    % ...and fourth column: Twomey-Markowski

tools.textbar([0, length(ii)]); % initialize progress textbar

%% iterate through test cases, performing inversion and plotting results %%

for i = ii

    % import non-inverted (a.k.a raw) SMPS files and average particle counts
    dat(i) = DAT(i, fname, tname, dm_rowid(i));
    
    % acquire and save the folder address for the TSI's SMPS data files 
    fadd_tsi_inv = intab.fdir_tsi_inv{i};
    ii0 = i;
    while isempty(fadd_tsi_inv) && (ii0 > 0)
        ii0 = ii0 - 1;
        fadd_tsi_inv = intab.fdir_tsi_inv{ii0};
    end
    fadd_tsi_noninv = intab.fdir_tsi_noninv{i};
    ii0 = i;
    while isempty(fadd_tsi_noninv) && (ii0 > 0)
        ii0 = ii0 - 1;
        fadd_tsi_noninv = intab.fdir_tsi_noninv{ii0};
    end

    % identify the address of TSI files
    fadd_tsi_inv = strcat(fadd_tsi_inv, '\',...
        intab.fname_tsi{i}, '.csv');
    fadd_tsi_noninv = strcat(fadd_tsi_noninv, '\',...
        intab.fname_tsi{i}, '.csv');
    
    % import the size distribution results of TSI software
    dist_tsi_inv(i) = hn.import_dist(fadd_tsi_inv, dat(i).sid,...
        dat(i).df, dm_rowid(i)); % multiple charge and diffusion loss corrected
    dist_tsi_noninv(i) = hn.import_dist(fadd_tsi_noninv, dat(i).sid,...
        dat(i).df, dm_rowid(i)); % only diffusion loss corrected
    
    % adjust the description text
    % rigstr{i} = DAT.LABELADJ(rigstr{i});
    
    % reconstructed mobility points
    dist_odias(i).d = logspace(log10(d_lim_1(i)), log10(d_lim_2(i)), 500)';

    % correct sheath flows
    prop.Q_c = dat(i).Qsh_dma;
    prop.Q_m = dat(i).Qsh_dma;

    % build inversion kernel
    A = kernel.gen_smps((dat(i).dm)', dist_odias(i).d, [], prop);
    
    % assign mean and std of particle counts
    b = (dat(i).dN)';
    Lb = bsxfun(@rdivide, eye(length(b)), (dat(i).sigma)');

    disp(' ');

    % Twomey-Markowski's inversion
    xi = invert.get_init(Lb * A, Lb * b, dist_odias(i).d, dat(i).dm);
    dist_odias(i).x_twomark = invert.twomark(Lb * A, Lb * b, length(xi),...
        xi);

    % correct Twomey-Markowski for concentration
    rN(i,4) = dist_tsi_inv(i).A_tot /...
        trapz(log10(dist_odias(i).d((dist_odias(i).d <= max(dat(i).dm))...
        & (dist_odias(i).d >= min(dat(i).dm)))),...
        dist_odias(i).x_twomark((dist_odias(i).d <= max(dat(i).dm)) &...
        (dist_odias(i).d >= min(dat(i).dm))));
    dist_odias(i).x_twomark = rN(i,4) * dist_odias(i).x_twomark;

    % get geometric standard deviation from Twomey-Markowski for later...
        % ...use in exponential distance
    w_twomark = dist_odias(i).x_twomark / sum(dist_odias(i).x_twomark);
    d_gm_twomark = 10^(sum(w_twomark .* log10(dist_odias(i).d))); % GM
    sigma_g_twomark = 10^(sqrt(sum(w_twomark .* (log10(dist_odias(i).d) -...
        log10(d_gm_twomark)).^2))); % GSD

    % perform inversion using 1st-order Tikhonov
    [dist_odias(i).x_tk1, ~, ~, Gpo_inv_tk1] = ...
        invert.tikhonov(Lb * A, Lb * b, lambda_tk1(i), 1);
    dist_odias(i).Gpo_tk1 = inv(Gpo_inv_tk1);
    
    % get and apply the correction factor for concentration in Tikhonov's...
        % ...1st-order method
    rN(i,1) =  dist_tsi_inv(i).A_tot /...
        trapz(log10(dist_odias(i).d((dist_odias(i).d <= max(dat(i).dm)) &...
        (dist_odias(i).d >= min(dat(i).dm)))),...
        dist_odias(i).x_tk1((dist_odias(i).d <= max(dat(i).dm)) &...
        (dist_odias(i).d >= min(dat(i).dm))));
    dist_odias(i).x_tk1 = rN(i,1) * dist_odias(i).x_tk1;
    dist_odias(i).Gpo_tk1 = rN(i,1) .* dist_odias(i).Gpo_tk1;

    % invert with 2nd-order Tikhonov
    [dist_odias(i).x_tk2, ~, ~, Gpo_inv_tk2] = ...
        invert.tikhonov(Lb * A, Lb * b, lambda_tk2(i), 2);
    dist_odias(i).Gpo_tk2 = inv(Gpo_inv_tk2);

    % correct 2nd-order Tikhonov for concentration
    rN(i,2) =  dist_tsi_inv(i).A_tot /...
        trapz(log10(dist_odias(i).d((dist_odias(i).d <= max(dat(i).dm)) &...
        (dist_odias(i).d >= min(dat(i).dm)))),...
        dist_odias(i).x_tk2((dist_odias(i).d <= max(dat(i).dm)) &...
        (dist_odias(i).d >= min(dat(i).dm))));
    dist_odias(i).x_tk2 = rN(i,2) * dist_odias(i).x_tk2;
    dist_odias(i).Gpo_tk2 = rN(i,2) .* dist_odias(i).Gpo_tk2;    

    % invert with exponential distance method
    [dist_odias(i).x_ed, ~, ~, Gpo_inv_ed] = ...
        invert.exp_dist(Lb * A, Lb * b, lambda_ed(i),...
        log10(sigma_g_twomark), dist_odias(i).d);
    dist_odias(i).Gpo_ed = inv(Gpo_inv_ed);

    % correct exponential distance results for concentration
    rN(i,3) =  dist_tsi_inv(i).A_tot /...
        trapz(log10(dist_odias(i).d((dist_odias(i).d <= max(dat(i).dm)) &...
        (dist_odias(i).d >= min(dat(i).dm)))),...
        dist_odias(i).x_ed((dist_odias(i).d <= max(dat(i).dm)) &...
        (dist_odias(i).d >= min(dat(i).dm))));
    dist_odias(i).x_ed = rN(i,3) * dist_odias(i).x_ed;
    dist_odias(i).Gpo_ed = rN(i,3) .* dist_odias(i).Gpo_ed;    

    if strcmp(opts2.correct, 'on') || strcmp(opts2.correct, 'On') ||...
            strcmp(opts2.correct, 'ON')

        % initialize a plot to check inversion with baseline data at...
        % ...every iteration
        f2 = figure(2);
        f2.Position = [200, 200, 550, 650];
        set(f2, 'color', 'white');        
        
        chk_flag1 = true; % flag for the loop that reiterates the inversion
        chk_flag2 = false;
        
        while chk_flag1
            
            figure(f2) % navigate to the plot that shows iterations on inversion

            if chk_flag2
                
                % rebuild kernel
                dist_odias(i).d = logspace(log10(d_lim_1(i)),...
                    log10(d_lim_2(i)), 500)';
                A = kernel.gen_smps((dat(i).dm)', dist_odias(i).d, [], prop);

                % repeat Twomey-Markowski (d might have changed)
                xi = invert.get_init(Lb * A, Lb * b, dist_odias(i).d,...
                    dat(i).dm);
                dist_odias(i).x_twomark = invert.twomark(Lb * A,...
                    Lb * b, length(xi), xi);
                rN(i,4) = dist_tsi_inv(i).A_tot /...
                    trapz(log10(dist_odias(i).d((dist_odias(i).d <=...
                    max(dat(i).dm)) & (dist_odias(i).d >= min(dat(i).dm)))),...
                    dist_odias(i).x_twomark((dist_odias(i).d <=...
                    max(dat(i).dm)) & (dist_odias(i).d >= min(dat(i).dm))));
                dist_odias(i).x_twomark = rN(i,4) *...
                    dist_odias(i).x_twomark;

                % repeat 1s-order Tikhonov with new smoothing factor
                [dist_odias(i).x_tk1, ~, ~, Gpo_inv_tk1] = ...
                    invert.tikhonov(Lb * A, Lb * b, lambda_tk1(i), 1);
                dist_odias(i).Gpo_tk1 = inv(Gpo_inv_tk1);
                rN(i,1) = dist_tsi_inv(i).A_tot /...
                    trapz(log10(dist_odias(i).d((dist_odias(i).d <=...
                    max(dat(i).dm)) & (dist_odias(i).d >= min(dat(i).dm)))),...
                    dist_odias(i).x_tk1((dist_odias(i).d <= max(dat(i).dm)) &...
                    (dist_odias(i).d >= min(dat(i).dm))));
                dist_odias(i).x_tk1 = rN(i,1) * dist_odias(i).x_tk1;
                dist_odias(i).Gpo_tk1 = rN(i,1) .* dist_odias(i).Gpo_tk1;

                % repeat 2nd-order Tikhonov with new smoothing factor
                [dist_odias(i).x_tk2, ~, ~, Gpo_inv_tk2] = ...
                    invert.tikhonov(Lb * A, Lb * b, lambda_tk2(i), 2);
                dist_odias(i).Gpo_tk2 = inv(Gpo_inv_tk2);
                rN(i,2) = dist_tsi_inv(i).A_tot /...
                    trapz(log10(dist_odias(i).d((dist_odias(i).d <=...
                    max(dat(i).dm)) & (dist_odias(i).d >= min(dat(i).dm)))),...
                    dist_odias(i).x_tk2((dist_odias(i).d <= max(dat(i).dm)) &...
                    (dist_odias(i).d >= min(dat(i).dm))));
                dist_odias(i).x_tk2 = rN(i,2) * dist_odias(i).x_tk2;
                dist_odias(i).Gpo_tk2 = rN(i,2) .* dist_odias(i).Gpo_tk2;

                % repeat exponential distacne with new smoothing factor
                [dist_odias(i).x_ed, ~, ~, Gpo_inv_ed] = ...
                    invert.exp_dist(Lb * A, Lb * b, lambda_ed(i),...
                    log10(sigma_g_twomark), dist_odias(i).d);
                dist_odias(i).Gpo_ed = inv(Gpo_inv_ed);
                rN(i,3) = dist_tsi_inv(i).A_tot /...
                    trapz(log10(dist_odias(i).d((dist_odias(i).d <=...
                    max(dat(i).dm)) & (dist_odias(i).d >= min(dat(i).dm)))),...
                    dist_odias(i).x_ed((dist_odias(i).d <= max(dat(i).dm)) &...
                    (dist_odias(i).d >= min(dat(i).dm))));
                dist_odias(i).x_ed = rN(i,3) * dist_odias(i).x_ed;
                dist_odias(i).Gpo_ed = rN(i,3) .* dist_odias(i).Gpo_ed;

            end
                        
            % plot Twomey-Markowski's output (which is closest to TSI's inversion)
            plt2_twomark = hn.plotci_custom(dist_odias(i).d,...
                dist_odias(i).x_twomark, [], [], hex2rgb(clr2{1}), '-');
            hold on
            
            % plot TSI's non-corrected data
            plt2_tsi_noninv = hn.plotci_custom(dist_tsi_noninv(i).dm',...
                dist_tsi_noninv(i).dn_dlogdm', [], [], hex2rgb(clr2{2}), '-.');
            hold on
            
            % plot Tikhonov's 1st-order results based on the last assigned lambda
            plt2_tk1 = hn.plotci_custom(dist_odias(i).d, dist_odias(i).x_tk1,...
                dist_odias(i).Gpo_tk1, [], hex2rgb(clr2{4}), '-');
            hold on

            % plot exponential distance's updated results
            plt2_ed = hn.plotci_custom(dist_odias(i).d, dist_odias(i).x_ed,...
                dist_odias(i).Gpo_ed, [], hex2rgb(clr2{6}), '-');
            hold on

            % plot TSI's multiple charge corrected data
            plt2_tsi_inv = hn.plotci_custom(dist_tsi_inv(i).dm',...
                dist_tsi_inv(i).dn_dlogdm', [], [], hex2rgb(clr2{3}), '-');
            hold on

            % plot Tikhonov's updated 2nd-order results
            plt2_tk2 = hn.plotci_custom(dist_odias(i).d, dist_odias(i).x_tk2,...
                dist_odias(i).Gpo_tk2, [], hex2rgb(clr2{5}), '-.');
            
            % set the appearance
            set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
                'TickLength', [0.02 0.02], 'XScale', 'log')
            xlim([d_lim_1(i), d_lim_2(i)])
            xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
                'FontSize', 16)
            ylabel('$\mathrm{d}n/\mathrm{dlog}(d_\mathrm{m} [\#/\mathrm{cm}^3])$',...
                'interpreter', 'latex', 'FontSize', 16)
            title(rigstr{i}, 'interpreter', 'latex', 'FontSize', 12)
            subtitle(' ', 'interpreter', 'latex', 'FontSize', 12)
            legend([plt2_tk1, plt2_tk2, plt2_ed, plt2_twomark,...
                plt2_tsi_inv, plt2_tsi_noninv],...
                {'Tikhonov 1st-order', 'Tikhonov 2nd-order',...
                'Exponential distance', 'Twomey-Markowski',...
                'TSI-MCC', 'TSI-nonMCC'}, 'interpreter', 'latex',...
                'FontSize', 12, 'Location', 'southoutside', 'orientation',...
                'horizontal', 'NumColumns', 2)
            
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
                        msg_params_refine = {strcat('Please enter a refined',...
                            string(newline), 'lambda_tk1:'), 'lambda_tk2:',...
                            'lambda_ed:', 'dm_min:', 'dm_max:'};
                        chk_flag2 = true; % flag to redo the inversion
                    else                        
                        msg_params_refine = {strcat('Please enter a refined',...
                            string(newline), 'lambda_tk1: (initial value:',...
                            {' '}, num2str(lambda_tk1_0(i), '%.1e'), ',',...
                            {' '}, 'recent value:', {' '},...
                            num2str(lambda_tk1(i), '%.1e'), ')'),...
                            cell2mat(strcat('lambda_tk2: (initial value:',...
                            {' '}, num2str(lambda_tk2_0(i), '%.1e'), ',',...
                            {' '}, 'recent value:', {' '},...
                            num2str(lambda_tk2(i), '%.1e'), ')')),...                            
                            cell2mat(strcat('lambda_ed: (initial value:',...
                            {' '}, num2str(lambda_ed_0(i), '%.1e'), ',',...
                            {' '}, 'recent value:', {' '},...
                            num2str(lambda_ed(i), '%.1e'), ')')),...
                            cell2mat(strcat('dm_min:', {' '}, '(initial value:',...
                            {' '}, num2str(d_lim_1_0(i), '%.1e'), ',', {' '},...
                            'recent value:', {' '}, num2str(d_lim_1(i), '%.1e'),...
                            ')')), cell2mat(strcat('dm_max:', {' '}, '(initial value:',...
                            {' '}, num2str(d_lim_2_0(i), '%.1e'), ',', {' '},...
                            'recent value:', {' '}, num2str(d_lim_2(i), '%.1e'),...
                            ')'))};
                    end
                    
                    % ask user to revise lambda
                    response2b = inputdlg(msg_params_refine,...
                        'Refine Tikhonov parameters', [1 60; 1 60; 1 60;...
                        1 60; 1 60], {num2str(lambda_tk1(i)),...
                        num2str(lambda_tk2(i)), num2str(lambda_ed(i)),...
                        num2str(d_lim_1(i)), num2str(d_lim_2(i))});
                    
                    % terminate the code if user hits cancel
                    if ~isempty(response2b)
                        lambda_tk1(i) = str2double(response2b{1});
                        lambda_tk2(i) = str2double(response2b{2});
                        lambda_ed(i) = str2double(response2b{3});
                        d_lim_1(i) = str2double(response2b{4});
                        d_lim_2(i) = str2double(response2b{5});
                    else
                        return;
                    end
                    
                case 'Cancel'
                    % terminate the program
                    return
            end
            
        end
        
        % extents of TSI mobility setpoints
        dmax_i = max(dat(i).dm);
        dmin_i = min(dat(i).dm);
    
        % indices of the constructed setpoints that fall within the...
            % ...experimental setpoints
        inds_in = (dist_odias(i).d <= dmax_i) & (dist_odias(i).d >= dmin_i);
    
        % determine final inversion to be output for within the range of TSI...
            % ...detection    
        response2c = listdlg('PromptString', {'Select method for output within',...
            'the TSI detection range', ' '}, 'SelectionMode','single',...
            'ListString', {'Tikhonov 1st order', 'Tikhonov 2nd order',...
            'Exponential distance'});
        
        % terminate the program if hit cancel
        if isempty(response2c); return; end
        
        switch response2c
            case 1
                dist_odias(i).x(inds_in) = dist_odias(i).x_tk1(inds_in);
                Gpo0 = dist_odias(i).Gpo_tk1;
                Gpo00 = dist_odias(i).Gpo_tk1(inds_in,inds_in);
                
            case 2
                dist_odias(i).x(inds_in) = dist_odias(i).x_tk2(inds_in);
                Gpo0 = dist_odias(i).Gpo_tk2;
                Gpo00 = dist_odias(i).Gpo_tk2(inds_in,inds_in);
    
            case 3
                dist_odias(i).x(inds_in) = dist_odias(i).x_ed(inds_in);
                Gpo0 = dist_odias(i).Gpo_ed;
                Gpo00 = dist_odias(i).Gpo_ed(inds_in,inds_in);
                
        end
    
        % indices of the constructed setpoints that fall outside the...
            % ...experimental setpoints    
        inds_out = (dist_odias(i).d > dmax_i) | (dist_odias(i).d < dmin_i);
            
        response2d = listdlg('PromptString', {'Select method for output outside',...
            'the TSI detection range', ' '}, 'SelectionMode','single',...
            'ListString', {'Tikhonov 1st order', 'Tikhonov 2nd order',...
            'Exponential distance', 'Twomey-Markowski'});
        
        % terminate the program if hit cancel
        if isempty(response2d); return; end
        
        switch response2d
            case 1
                dist_odias(i).x(inds_out) = dist_odias(i).x_tk1(inds_out);
                dist_odias(i).Gpo = dist_odias(i).Gpo_tk1;
                
            case 2
                dist_odias(i).x(inds_out) = dist_odias(i).x_tk2(inds_out);
                dist_odias(i).Gpo = dist_odias(i).Gpo_tk2;
                
            case 3
                dist_odias(i).x(inds_out) = dist_odias(i).x_ed(inds_out);
                dist_odias(i).Gpo = dist_odias(i).Gpo_ed;
                
            case 4
                dist_odias(i).x(inds_out) = dist_odias(i).x_twomark(inds_out);
                
        end
    
        % terminate the program if hit cancel
        if isempty(response2d); return; end    
        
        if response2d == 4
            dist_odias(i).Gpo = Gpo0;
        
        else
            dist_odias(i).Gpo(inds_in,inds_in) = Gpo00;
        end
    
        % adjust the direction of output counts vector
        dist_odias(i).x = dist_odias(i).x'; 
        
        % clear figure for selecting inversion method
        clf(f2, 'reset')

    else

        % if not asking from user, use Tikhonov's 2nd order as final output
        dist_odias(i).x = dist_odias(i).x_tk2;
        dist_odias(i).Gpo = dist_odias(i).Gpo_tk2;
        
    end

    figure(f1) % navigate to final plot for size distribution
    nexttile(1) % non-log scales    
    
    % display the inverted mobility size distribution based on Tikhonov
    plt1{i} =  hn.plotci_custom(dist_odias(i).d, dist_odias(i).x,...
        dist_odias(i).Gpo, [], hex2rgb(clr1{i}), linstl{i});
    hold on
    
    % rigstr{i} = strrep(rigstr{i},'nof','n/f'); % minor label adjustment

    % apply the axis configurations
    if i == ii(end)
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
            'TickLength', [0.02 0.02], 'XScale', 'log')    
        xlim([15,950])
        xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
            'FontSize', 16)
        ylabel('$\mathrm{d}n/\mathrm{dlog}(d_\mathrm{m} [\#/\mathrm{cm}^3])$',...
            'interpreter', 'latex', 'FontSize', 16)
    end
    
    nexttile(2) % log scales
    
    % remove zero counts for plotting in log scale
    dm2{i} = dist_odias(i).d(dist_odias(i).x > 0);
    dN2{i} = dist_odias(i).x(dist_odias(i).x > 0);
    
    % plot nonzero scales for display in log-log scale
    opts1b.log = 'on';
    hn.plotci_custom(dm2{i}, dN2{i},...
        dist_odias(i).Gpo(dist_odias(i).x > 0, dist_odias(i).x > 0),...
        [], hex2rgb(clr1{i}), linstl{i}, opts1b)
    hold on

    % apply the log-log scale and other axis configurations
    if i == ii(end)
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
            'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')    
        xlim([15,950])
        ylim([-Inf 1.2 * max(cat(1,dN2{:}))])
        xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
            'FontSize', 16)
        ylabel('$\mathrm{d}n/\mathrm{dlog}(d_\mathrm{m}) [\#/\mathrm{cm}^3]$',...
            'interpreter', 'latex', 'FontSize', 16)

        % ask user to adjust y axis in log scale
        ymin_log = str2double(inputdlg('Please enter minimum of vertical axis in Log scale',...
            'Adjust axis', [1 50], {'1e3'}));
        if isempty(ymin_log)
            return % stop code execution if hit cancel
        else
            ylim([ymin_log 1.2 * max(cat(1,dN2{:}))]) % adjust y axis in...
                % ...log scale
        end
    end

    % find local modes of distribution and remove noise
    [dn_dlogd_mode, dist_odias(i).d_mode] = findpeaks(dist_odias(i).x_tk1,...
        dist_odias(i).d); 
    dist_odias(i).d_mode(dn_dlogd_mode / max(dn_dlogd_mode) < 0.1) = [];

    % find geometric mean (GM) and geometric standard deviation (GSD)
    w = dist_odias(i).x_tk1 / sum(dist_odias(i).x_tk1); % normalize dn/dlog(d) to get weights
    dist_odias(i).d_gm = 10^(sum(w .* log10(dist_odias(i).d))); % GM
    dist_odias(i).sigma_g = 10^(sqrt(sum(w .* (log10(dist_odias(i).d) -...
        log10(dist_odias(i).d_gm)).^2))); % GSD

    % total concentration (i.e. area below the size distribution curve)
    dist_odias(i).n_tot = trapz(log10(dist_odias(i).d), dist_odias(i).x_tk1);
    
    % skewness in log-space (0 for a normal distribution, positive for...
    % ...right-skewed, negative for left-skewed)
    dist_odias(i).skw = sum(w .* (log10(dist_odias(i).d) -...
        log10(dist_odias(i).d_gm)).^3) / ((log10(dist_odias(i).sigma_g))^3);
    
    % skewness in log-space (3 for a normal distribution, > 3 for...
    % ...heavy tails, < 3 for light tails)
    dist_odias(i).krts = sum(w .* (log10(dist_odias(i).d) -...
        log10(dist_odias(i).d_gm)).^4) / ((log10(dist_odias(i).sigma_g))^4);

    tools.textbar([find(ii == i), length(ii)]); % update textbar
end

lgd1 = legend(cat(2, plt1{:}), cat(2, rigstr(ii)), 'interpreter', 'latex',...
    'FontSize', 12);
lgd1.Layout.Tile = 'south';
lgd1.NumColumns = 2;

% close the inversion iteration figure if present
if exist('f2', 'var'); close(f2); end

%% present data in groups if requested by user

if exist('grp', 'var') && ~isempty(grp)
    
    % initialize variables to process and store averages
    d0 = cell(n_grp,1);
    dn0_dlogd = cell(n_grp,1);
    d = cell(n_grp,1);
    dn_dlogd = cell(n_grp,1);
    se_dn_dlogd = cell(n_grp,1);
    ub_dn_dlogd = cell(n_grp,1);
    lb_dn_dlogd = cell(n_grp,1);

    % initialize size distribution figure for groups
    f3 = figure(3);
    f3.Position = [150, 150, 500, 500];
    set(f3, 'color', 'white');
    
    % plot and legend placeholders
    plt3 = cell(n_grp,1);
    lgdtxt3 = cell(n_grp,1);

    % loop over the groups
    for j = 1 : n_grp

        d0{j} = cell(length(ind_grp{j}),1);
        dn0_dlogd{j} = cell(length(ind_grp{j}),1);

        for k = 1 : length(ind_grp{j})

            d0{j}{k} = dist_odias(ind_grp{j}(k)).d;
            dn0_dlogd{j}{k} = dist_odias(ind_grp{j}(k)).x;

        end
        
        % take average of size distribution over each group
        [dn_dlogd{j}, d{j}, se_dn_dlogd{j}] = hn.avg_dist(d0{j},...
            dn0_dlogd{j});

        % compute the bounds
        ub_dn_dlogd{j} = dn_dlogd{j} + se_dn_dlogd{j};
        lb_dn_dlogd{j} = dn_dlogd{j} - se_dn_dlogd{j};
        
        % plot the average distributions and error bounds
        plt3{j} = plot(d{j}, dn_dlogd{j}, 'Color', hex2rgb(clr3{j}{1}),...
            'LineWidth', 1.5);
        hold on
        fill([d{j}, fliplr(d{j})], [ub_dn_dlogd{j}, fliplr(lb_dn_dlogd{j})],...
             hex2rgb(clr3{j}{2}), 'EdgeColor', 'none', 'FaceAlpha', 0.5);

        lgdtxt3{j} = strcat('Group', {' '}, num2str(j));

    end
    
    % set apprearances
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
        'TickLength', [0.02 0.02], 'XScale', 'log')
    xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
        'FontSize', 16)
    ylabel('$\mathrm{d}n/\mathrm{dlog}(d_\mathrm{m}) [\#/\mathrm{cm}^3]$',...
        'interpreter', 'latex', 'FontSize', 16)
    xlim([15,1000])
    legend(cat(2, plt3{:}), cat(2, lgdtxt3{:}), 'interpreter', 'latex',...
    'FontSize', 12, 'Location', 'northwest');
    
end

%% save figures and data

% make directory for results to be saved
fdir_out = strcat('outputs\', tname0);
if ~isfolder(fdir_out); mkdir(fdir_out); end

% customize the name of directory
fname_out = strcat(tname0, '_', datestr(datetime('now'))); %#ok<DATST>
fname_out = regexprep(fname_out, ':', '-');
fname_out = regexprep(fname_out, ' ', '_');

% save individual inverted size distributions
saveas(f1, strcat(fdir_out, '\Fig1_', fname_out, '.fig'));

% save average size distributions
if exist('f3', 'var')
    saveas(f3, strcat(fdir_out, '\Fig3_', fname_out, '.fig'));
end

% save MATLAB worspace
save(strcat(fdir_out, '\', fname_out, '.mat'));


