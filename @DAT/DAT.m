classdef DAT
    % a class of mobility size distribution data obtained from DMA scan...
    %   ...files (typically used per AAC setpoints in tandem...
    %   ...measurements but also capable of handling SMPS-only scans)
    
    properties

        f_dir = [] % directory to the DMA data file
        f_name = [] % name of the DMA data file
        da = [] % AAC setpoint correspodning to the DMA scan
        Qsh_aac = [] % AAC sheath flow rate per AAC setpoint 
        Qsh_dma = [] % DMA sheath flow rate per AAC setpoint 
        dm = [] % DMA scan setpoints
        dm_max = [] % mode of mobility size distribution
        n_dm = [] % number of DMA setpoints per AAC setpoint
        dN0 = [] % raw DMA counts
        sid = [] % DMA scan ids to be extracted
        df = [] % dilution factor
        dN = [] % averaged DMA counts
        sigma = [] % a vector of sd of counts in each bin
        lbl = [] % strcuture containing data labels for visualization

        % the following are only used as placeholders of original data...
        %   ...while using the COMBINE function
        Qsh0_aac = []
        Qsh0_dma = []
        dm0 = []
        n0_dm = []

    end
    
    methods
        
        function obj = DAT(pid, f0_name, t_name, dm_rowid, varargin)
            % DAT constructs an instance of this class.
            %   DMA .xlsx file is imported. Raw counts and setpoint info...
            %   ...are extracted.
            % ======================================================== %
            % f0_name: filename of instrument settings table
            % t_name: tab (sheet) name in excel file
            % dm_rowid: assign which row of DMA file to use for...
            %   ...extracting mobility diameter.
            % pid: the test point id for which DMA data is to be extracted 
            % varargin: variable argument; varargin{1} is address of...
            %   ...instrument settings file.
            % ======================================================== %
            
            if nargin == 0; return; end % just initialize class

            % import table of instrument settings
            if nargin < 4; varargin = []; end
            if (isempty(varargin)); varargin = 'inputs'; end
            f0_add = char(strcat(varargin, {'\'}, f0_name, '.xlsx'));
            instab = readtable(f0_add, 'UseExcel', true, 'Sheet', t_name);
            
            % find and save the file directory in the table
            obj.f_dir = instab.f_add_raw{pid};
            pid0 = pid;
            while isempty(obj.f_dir) && (pid0 > 0)
                pid0 = pid0 - 1;
                obj.f_dir = instab.f_add_raw{pid0};
            end
            
            % save the file name, AAC setpoint, AAC and DMA sheath flows...
                % ...and DMA scan indices
            obj.f_name = instab.f_name_raw{pid};
            flds = fieldnames(instab);
            if ismember({'da_nm_'},flds)
                obj.da = instab.da_nm_(pid);
            end
            if ismember({'Qsh_AAC_lpm_'},flds)
                obj.Qsh_aac = instab.Qsh_AAC_lpm_(pid) / 60000; % convert to m3/s
            end
            obj.Qsh_dma = instab.Qsh_DMA_lpm_(pid) / 60000; % convert to m3/s
            
            % load the DMA data
            f_add = char(strcat(obj.f_dir, {'\'}, obj.f_name, '.csv'));
            opts = detectImportOptions(f_add);
            opts = setvartype(opts, 'double');
            dat0 = readtable(f_add, opts);
            
            % save mobility setpoints and their number
            obj.dm = table2array(dat0(dm_rowid, 39:end));
            obj.n_dm = length(obj.dm);

            % save raw counts for the selected scans
            dN_raw = table2array(dat0((dm_rowid+1):end, 39:end));
            sidfun = str2func(['@(x)', instab.sid{pid}]);
            obj.sid = sidfun(); % scan indices to be selected
            obj.df = instab.DF(pid); % dilition factor to be multiplied...
            %   ...by raw counts
            if isempty(obj.sid)
                obj.dN0 = dN_raw;
                obj.sid = 1 : size(dN_raw,1);
            else
                obj.dN0 = dN_raw(obj.sid,:); % extract selected scans                    
            end
            obj.dN0 = obj.df * obj.dN0; % correct counts for dilution

            % average and standard deviation of counts in each bin
            obj.dN = mean(obj.dN0);
            obj.sigma = max(std(obj.dN0), 1e-3 * max(std(obj.dN0)));
            
            % calculate the mode of a mobility size distribution
            obj.dm_max = obj.dm(obj.dN == max(obj.dN));
            
            if ~isempty(obj.da)
                % make a label for the object (using AAC setpoint)
                obj.lbl = {strcat('$d_\mathrm{a}$ = ',...
                    num2str(obj.da, '%.1f'), ' [nm]')};
            end

        end
        
        function obj = SELECT(obj, sid)
            % SELECT picks data for a selected set of DMA scan repeats,...
            %   ...averages in each bin between the repeats, and...
            %   ...calculates standard deviation of repeats for later...
            %   ...inversion analyses.
            % ======================================================== %
            % sid: a number vector of scans to be averaged
            % ======================================================== %
            
            obj.dN = mean(obj.dN0(sid,:));
            obj.sigma = max(std(obj.dN0(sid,:)),...
                1e-3 * max(std(obj.dN0(sid,:))));
            
        end

        function obj = COMBINE(obj1, obj2)
            % COMBINE merges raw DMA scan data with (possibly) different...
            %   ...mobility ranges/resolutions/etc. for later...
            %   ...post-processing. The objects input to "COMBINE" are...
            %   ...presumably produced from different scan files.
            % ======================================================== %
            % obj1,2: the DMA data class objects to be combined 
            % ======================================================== %
            
            obj = DAT(); % initialize the output object
            
            % check if the AAC setpoint is the same (technically should be)
            if (obj1.da == obj2.da)
                obj.da = obj1.da;
            else
                error("Different AAC setpoints are being merged!")
            end
            
            % save file addresses
            if iscell(obj1.f_name)
                obj.f_dir = [obj1.f_dir, {obj2.f_dir}];
                obj.f_name = [obj1.f_name, {obj2.f_name}];
                n1_sid = length(horzcat(obj1.sid{:}));
            else
                obj.f_dir = {obj1.f_dir, obj2.f_dir};
                obj.f_name = {obj1.f_name, obj2.f_name};
                n1_sid = length(obj1.sid);
            end

            % save sheath flows
            n2_sid = length(obj2.sid);
            if ~isempty(obj1.Qsh_aac)
                if length(obj1.Qsh0_aac) > 1
                    obj.Qsh0_aac = [obj1.Qsh0_aac, obj2.Qsh_aac];
                else
                    obj.Qsh0_aac = [obj1.Qsh_aac, obj2.Qsh_aac];
                end
                obj.Qsh_aac = (n1_sid * obj1.n_dm * obj1.Qsh_aac +...
                    n2_sid * obj2.n_dm * obj2.Qsh_aac) /...
                    (n1_sid * obj1.n_dm + n2_sid * obj2.n_dm);
            end
            if length(obj1.Qsh0_dma) > 1
                obj.Qsh0_dma = [obj1.Qsh0_dma, obj2.Qsh_dma];
            else
                obj.Qsh0_dma = [obj1.Qsh_dma, obj2.Qsh_dma];
            end
            obj.Qsh_dma = (n1_sid * obj1.n_dm * obj1.Qsh_dma +...
                n2_sid * obj2.n_dm * obj2.Qsh_dma) /...
                (n1_sid * obj1.n_dm + n2_sid * obj2.n_dm);
            
            % categorize + compile input mobility setpoints and raw counts
            [dm_a, I_a] = setdiff(obj1.dm, obj2.dm, 'stable');
            [dm_c, I_c] = setdiff(obj2.dm, obj1.dm, 'stable');
            [dm_b, I_b1, I_b2] = intersect(obj1.dm, obj2.dm);
            obj.dm0 = {dm_a, dm_b, dm_c};
            if iscell(obj1.dN0)
                dN0_b = vertcat(repmat(obj1.dN(:, I_b1), n1_sid, 1),...
                    obj2.dN0(:, I_b2));
                obj.dN0 = [obj1.dN0, {obj2.dN0}];
                obj.sid = [obj1.sid, {obj2.sid}];
                obj.n0_dm = [obj1.n0_dm, obj2.n_dm];
            else
                dN0_b = vertcat(obj1.dN0(:, I_b1), obj2.dN0(:, I_b2));
                obj.dN0 = {obj1.dN0, obj2.dN0};
                % obj.dN0 = {obj1.dN0(:, I_a), dN0_b, obj2.dN0(:, I_c)};
                obj.sid = {obj1.sid, obj2.sid};
                obj.n0_dm = [obj1.n_dm, obj2.n_dm];
            end
            obj.df = [obj1.df, obj2.df];

            % organize mobility setpoints and averages and reaverage...
            %   ...counts when needed
            obj.dm = horzcat(dm_a, dm_b, dm_c);
            [obj.dm, I] = sort(obj.dm);
            obj.n_dm = length(obj.dm);
            dN_b = mean(dN0_b);
            obj.dN = [obj1.dN(:, I_a), dN_b, obj2.dN(:, I_c)];
            obj.dN = obj.dN(I);
            if iscell(obj1.dN0)
                sigma_b = sqrt(((n1_sid - 1) * obj1.sigma(:, I_b1).^2 +...
                    (n2_sid - 1) * obj2.sigma(:, I_b2).^2 +...
                    n1_sid * (obj1.dN(:, I_b1) - dN_b).^2 +...
                    n2_sid * (obj2.dN(:, I_b2) - dN_b).^2) /...
                    (n1_sid + n2_sid - 1));
            else
                sigma_b = std(dN0_b);
            end
            sigma_b = max(sigma_b, 1e-3 * max(sigma_b));
            obj.sigma = [obj1.sigma(:, I_a), sigma_b, obj2.sigma(:, I_c)];
            obj.sigma = obj.sigma(I);

            obj.dm_max = obj.dm(obj.dN == max(obj.dN)); % update mode

            obj.lbl = obj1.lbl; % recall label

        end
        
        function objs = PAIRFINDER(objs)
            % PAIRFINDER looks for data objects with similar AAC...
            %   ...setpoints and combines them.
            % ======================================================== %
            % objs: the DMA data class objects to be merged 
            % ======================================================== %

            da0 = cat(2, objs.da);
            i0_pair = nchoosek(1:length(da0), 2);
            da0 = da0(i0_pair);
            i_pair = find(da0(:,1) == da0(:,2), 1);

            while ~isempty(i_pair)

                objs(i0_pair(i_pair,1)) = COMBINE(objs(i0_pair(i_pair,1)),...
                    objs(i0_pair(i_pair,2)));

                objs(i0_pair(i_pair,2)) = [];

                da0 = cat(2, objs.da);
                i0_pair = nchoosek(1:length(da0), 2);
                da0 = da0(i0_pair);
                i_pair = find(da0(:,1) == da0(:,2), 1);                
                
            end
            
        end

        function obj = COMPILE(objs)
            % COMPILE puts together data from a series of objects of...
            %   ...class DAT into a single output object for...
            %   ...post-processing purposes.
            % ======================================================== %
            % objs: the DMA data class objects to be merged 
            % ======================================================== %
            
            obj = DAT(); % initialize the output object
            
            % concatenate mobility-resolved properties
            obj.dm = cat(2, objs.dm);
            obj.dN = cat(2, objs.dN);
            obj.sigma = cat(2, objs.sigma);

            % concatenate non-mobility-resolved properties and match...
                % ...them with mobility-resolved ones
            obj.n_dm = cat(2, objs.n_dm);
            obj.da = cat(2, objs.da);
            if ~isempty(obj.da); obj.da = repelem(obj.da, obj.n_dm); end
            obj.Qsh_aac = cat(2, objs.Qsh_aac);
            if ~isempty(obj.Qsh_aac); obj.Qsh_aac = repelem(obj.Qsh_aac, obj.n_dm); end
            obj.Qsh_dma = cat(2, objs.Qsh_dma);
            obj.Qsh_dma = repelem(obj.Qsh_dma, obj.n_dm);
            obj.n_dm = repelem(obj.n_dm, obj.n_dm);            
            
        end

        function h = PLOTDIST_TAND(objs, varargin)
            % PLOTDIST_TAND displays a set of mobility size...
            %   ...distributions existing within the input object. It...
            %   ...also colorcodes them based on selected setpoints.
            % ======================================================== %
            % objs: a class of DAT objects
            % varargin: variable argument; varargin{1} is plot title.
            % h: output figure handle
            % ======================================================== %
            
            % initialize figure properties
            figure;
            h = gcf;
            h.Position = [0, 0, 900, 500];
            set(h, 'color', 'white');
            
            n = length(objs); % number of AAC setpoints

            % initialize plot colors
            cm = colormap(turbo); % choose colormap
            ii = round(1 + (length(cm) - 1) .* (0 : 1 / (n - 1) : 1));
            cm = cm(ii,:); % descretize colors
            % cm = flip(cm,1); % flip the order of colors if needed

            plt = cell(1,n); % placeholder for individual plot lines
            
            % plot the size distributions per AAC setpoint
            for i = 1 : n
                plt{i} = plot(objs(i).dm, objs(i).dN, 'Color', cm(i,:),...
                    'LineStyle', '-', 'LineWidth', 1);
                hold on
            end
            
            box on
            set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16,...
                'TickLength', [0.02 0.02], 'XScale', 'log',...
                'YScale', 'log')

            xlim([20, 1100])
            ylim([0.1, 1.1 * max(ylim)])

            xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
                'FontSize', 20)            
            ylabel('$\mathrm{d}N/Q_\mathrm{DMA}$ [$\mathrm{\#/cm\hat{A}^3}$]',...
                'interpreter', 'latex', 'FontSize', 20)

            legend(cat(2, plt{:}), cat(2, objs.lbl), 'interpreter',...
                'latex', 'FontSize', 16, 'Location', 'eastoutside',...
                'NumColumns', 2);
            
            if nargin > 1 && ~isempty(varargin{1})
                title(varargin{1}, 'FontSize', 16, 'interpreter','latex')
            end

        end
        
        function h = PLOTDIST(objs, varargin)
            % PLOTDIST displays and compares mobility size distributions...
            %   ...betweeen different experimental condistions.
            % ======================================================== %
            % objs: a class of DAT objects
            % varargin: variable argument; varargin{1} is plot title.
            % h: output figure handle
            % ======================================================== %
            
            % initialize variable input argument
            if nargin < 2
                if isempty(varargin) || ~isstruct(varargin{1})
                    varargin{1} = struct();
                end
            end

            % initialize titles describing experiments
            if ~isfield(varargin{1}, 'ttl')
                varargin{1}.ttl = '';
            end

            % initialize figure properties
            figure;
            h = gcf;
            h.Position = [0, 0, 900, 500];
            set(h, 'color', 'white');
            
            n = length(objs); % number of AAC setpoints

            % initialize plot colors
            cm = colormap(turbo); % choose colormap
            ii = round(1 + (length(cm) - 1) .* (0.1 : 0.8 / (n - 1) : 0.9));
            cm = cm(ii,:); % descretize colors
            % cm = flip(cm,1); % flip the order of colors if needed

            plt = cell(1,n); % placeholder for individual plot lines
            
            for i = 1 : n
                % plot average size distributions
                plt{i} = plot(objs(i).dm, objs(i).dN,...
                    'Color', cm(i,:), 'LineStyle', '-', 'LineWidth', 2.5);
                hold on

                % plot 95% confidence intervals
                plot(objs(i).dm, objs(i).dN + 1.96 * objs(i).sigma /...
                    sqrt(length(objs(i).sid)),...
                    objs(i).dm, objs(i).dN - 1.96 * objs(i).sigma /...
                    sqrt(length(objs(i).sid)),...
                    'Color', cm(i,:), 'LineStyle', ':', 'LineWidth', 1.5);
            end
            
            xlim([20, 1000])
            ylim([1, 1.5 * max(ylim)])
            
            box on
            set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16,...
                'TickLength', [0.02 0.02], 'XScale', 'log',...
                'YScale', 'log')

            xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
                'FontSize', 20)            
            ylabel('$\mathrm{d}N/Q_\mathrm{DMA}$ [$\mathrm{\#/cm\hat{A}^3}$]',...
                'interpreter', 'latex', 'FontSize', 20)
            
            if nargin > 1 && ~isempty(varargin{1}.ttl)
                legend(cat(2, plt{:}), cat(2, varargin{1}.ttl{:}),...
                    'interpreter', 'latex', 'FontSize', 16,...
                    'Location', 'eastoutside');            
            end

        end
        
        function h = DENSCATTER(objs, varargin)
            % DENSCATTER visualizes the individual counts in DMA bins...
            %   ...per AAC setpoints as a heatmap of scattered data.
            % ======================================================== %
            % objs: a class of DAT objects
            % varargin: variable argument; varargin{1} is the plotting...
            %   ...options variable.
            % h: output figure handle
            % ======================================================== %
            
            % initialize figure properties
            figure;
            h = gcf;
            h.Position = [0, 0, 500, 600];
            set(h, 'color', 'white');
            
            % define title, marker size and marker type if not input

            if nargin < 2
                if isempty(varargin) || ~isstruct(varargin{1})
                    varargin{1} = struct();
                end
            end

            % a cell array of titles describing each experiment
            if ~isfield(varargin{1}, 'ttl')
                varargin{1}.ttl = '';
            end

            % marker size
            if ~isfield(varargin{1}, 'ms')
                varargin{1}.ms = 80;
            end

            % marker type
            if ~isfield(varargin{1}, 'mt')
                varargin{1}.mt = 's';
            end

            % marker color
            mc = colormap(hot); % choose colormap
            % mc = flip(mc,1); % flip the order of colors if needed
                        
            % compile the properties of mobility bins over all objects
            obj = COMPILE(objs);

            % assign concentrations to colormap
            cdN = round(rescale(log(1 + obj.dN), 1, 256));
            cdN(isnan(cdN)) = 1;

            % plot the scattered heatmaps of mobility vs. aerodynamic...
                % ...diameter colored based on aerosol concentration 
            scatter(obj.dm, obj.da, varargin{1}.ms,...
                mc(cdN,:), varargin{1}.mt, 'filled');
            
            box on
            set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16,...
                'TickLength', [0.02 0.02], 'XScale', 'log',...
                'YScale', 'log')

            xlim([25, 1.05 * max(xlim)])
            ylim([38, 1.05 * max(ylim)])
            
            % print colorbar
            cb = colorbar;
            cb.Label.String = '$\mathrm{Log}\left(1 + \mathrm{d}N/\mathrm{d}N_\mathrm{max}\right)$';
            cb.Label.Interpreter  = 'latex';
            % cb.Label.Rotation = 360;
            cb.TickLabelInterpreter  = 'latex';
            cb.FontSize = 14;
            cb.Label.FontSize = 16;
            % cb.Limits = [1.0 1.6];
            % cb.Ticks = 1.0:0.1:1.6;
            cb.TickLength = 0.02;
            cb.Location = 'southoutside';
            % cbpos = get(cb, 'Position');
            % cb.Label.Position = [cbpos(1) - 0.75 , cbpos(2) + 1.705];
            
            xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
                'FontSize', 20)            
            ylabel('$$d_\mathrm{a}$ [nm]', 'interpreter',...
                'latex', 'FontSize', 20)
            
            if nargin > 1 && ~isempty(varargin{1}.ttl)
                title(varargin{1}.ttl, 'FontSize', 16,...
                    'interpreter','latex')
            end

        end
        
    end
    
    methods(Static)

        function cc_d = CC(d)
        
            % CC calculates a correction factor for Stokes' particle drag...
            %   ...equation based on Cunningham's suggestion.
            % =================================================================== %
            % d: characteristic diameter of particle
            % cc_d: correction factor
            % =================================================================== %
        
            % Cunningham's relation constants
            alpha = 1.165;
            beta = 0.483;
            gamma = 0.997;
            
            % calculate kinetic Knusden number
            lambda = 68; % mean free path of air @ room temperature [nm]
            kn_d = 2 * lambda ./ d;
            
            % calculate Cunningham's correction factor
            cc_d = 1 +  kn_d .* (alpha + beta * exp(-gamma ./ kn_d));
        
        end

        function rho_eff = RHO_EFF(dm, da)
            % RHO_EFF calculates effective density of a particle using...
            %   ...its aerodynamic and mobility diameters.
            % ======================================================== %
            % dm: mobility diameter
            % da: aerodynamic diameter
            % ======================================================== %
            
            % Cunningham's correction factors
            cc_da = DAT.CC(da);
            cc_dm = DAT.CC(dm);
            
            % calculate effective density
            rho_0 = 1000; % [kg/m3]; reference density
            rho_eff = rho_0 * (cc_da ./ cc_dm) .* (da ./ dm).^2;
            
        end

        function h = PLOTDENS(rho_eff, dm, varargin)
            % PLOTDENS plots a comparative form of 1d effective density...
            %   ...data vs. mobility diameter. It also colorcodes them...
            %   ...based on selected setpoints.
            % ======================================================== %
            % rho_eff: a cell array of effective densities for different...
            %   ...experiments
            % dm: corresponding modes of mobility diameter 
            % varargin: variable argument; varargin{1} is the plotting...
            %   ...options variable.            
            % h: output figure handle
            % ======================================================== %
            
            % initialize figure properties
            figure;
            h = gcf;
            h.Position = [0, 0, 900, 500];
            set(h, 'color', 'white');
            
            % get the number of experiments
            m = length(rho_eff);
            
            % initialize plotting options if not input
            if nargin > 2
                if isempty(varargin{1}) || ~isstruct(varargin{1})
                    varargin{1} = struct();
                end
                
                % a cell array of titles describing each experiment
                if ~isfield(varargin{1}, 'ttl')
                    varargin{1}.ttl = cell(1,m);
                end
                    
                % marker size
                if ~isfield(varargin{1}, 'ms')
                    varargin{1}.ms = 15 * ones(1,m);
                end

                % marker type 
                if ~isfield(varargin{1}, 'mt')
                    mt0 = {'^', 'o', 'v', 's', 'd', 'p' , '>', '<',...
                        'h', 'x', '+', '*', '_', '|'};
                    if m < 15
                        varargin{1}.mt = mt0(1:m);
                    else
                        error('Marker types need to be manually input!')
                    end
                end

            end

            % set marker colors
            mc = colormap(turbo); % choose colormap
            ii = round(1 + (length(mc) - 1) .* (0.1 : 0.8 / (m - 1) : 0.9));
            mc = mc(ii,:); % descretize colors
            mc = flip(mc,1); % flip the order of colors if needed
                        
            % placeholder for datasets corresponding to each experiments
            plt = cell(1, m+1);
            
            % plot the universal correlation of Olfert & Rogak (2019)
            r0 = (1e4 / 1e1)^(1 / (1e4 - 1));
            dm0 = 1e1 * ones(1e4, 1);
            for i = 2 : length(dm0)
                dm0(i) = dm0(i) * r0^(i-1);
            end
            rho0 = 510 * (dm0 / 100).^(-0.52);
            plt{m+1} = plot(dm0, rho0, 'Color', [0.4940 0.1840 0.5560],...
                'LineStyle', '-.', 'LineWidth', 3);
            hold on
            
            % plot the size distributions per AAC setpoint
            for i = 1 : m
                plt{i} = scatter(dm{i}, rho_eff{i}, varargin{1}.ms(i),...
                    mc(i,:), varargin{1}.mt{i});
            end
            
            box on
            set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16,...
                'TickLength', [0.02 0.02], 'XScale', 'log',...
                'YScale', 'log')

            xlim([50, 1000])
            ylim([100, 1200])

            xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
                'FontSize', 20)            
            ylabel('$\rho_\mathrm{eff}$ [$\mathrm{kg/m^3}$]', 'interpreter',...
                'latex', 'FontSize', 20)

            legend(cat(2, plt{:}), cat(2, varargin{1}.ttl{:},...
                {'Olfert and Rogak (2019)'}), 'interpreter', 'latex',...
                'FontSize', 14, 'Location', 'eastoutside');
            
        end
        
        function rigstr = LABELMAKER(f_name, t_name, varargin)
            % LABELMAKER makes text strings for use in data visualization.
            % ======================================================== %
            % f_name: filename for the table of properties describing...
            %   ...the experiment.
            % t_name: tab (sheet) name to be read within the excel file
            % varargin: variable argument; varargin{1} is directory...
            %   ...address of input table of properties.
            % rigstr: a 'string' variable that can be input to DAT...
            %   ...objects to generate data labels for plotting.
            % ======================================================== %
            
            % import table of test rig settings
            if nargin < 3; varargin = []; end
            if isempty(varargin); varargin = 'inputs'; end
            f_add = char(strcat(varargin, {'\'}, f_name, '.xlsx'));
            opts = detectImportOptions(f_add);
            opts = setvartype(opts, 'char');
            rigtab = readtable(f_add, opts, 'UseExcel',true,...
                'Sheet', t_name);
            
            % remove empty rows
            for i = 1 : size(rigtab,1)
                if isempty(cat(2,cell2mat(rigtab{i,:})))
                    i_rmv = i;
                    break
                end
            end

            if exist("i_rmv", "var"); rigtab(i_rmv:end, :) = []; end

            % extract field names from the table
            fld = fieldnames(rigtab);
            fld = fld(1:end-3);

            rigstr = cell(size(rigtab,1) - 1,1); % initialize output text
            
            for i = 1 : (size(rigtab,1) - 1)

                rigstr{i} = '';

                for j = 1 : length(fld)
    
                    if j > 1; rigstr{i} = strcat(rigstr{i}, {' '}); end
    
                    % print field name and value depending on being text or number
                    if isempty(str2num(cell2mat(rigtab{i,j}))) %#ok<ST2NM>
                        rigstr{i} = strcat(rigstr{i}, fld{j}, {': '}, rigtab{i,j});
                    else
                        rigstr{i} = strcat(rigstr{i}, {'$'}, fld{j}, {'$'}, {' = '},...
                            rigtab{i,j});
                    end
    
                    % print units
                    if ~isempty(cell2mat(rigtab{end,j}))
                        rigstr{i} = strcat(rigstr{i}, {' ['}, rigtab{end,j}, {']'});
                    end
    
                    % separate from next input
                    if j < length(fld)
                        rigstr{i} = strcat(rigstr{i}, {','});
                    end
    
                    rigstr{i} = cell2mat(rigstr{i});
    
                end
            end
            
        end

        function rigstr = LABELADJ(rigstr0, varargin)
            % LABELADJ adjusts the length and format of input label.
            % ======================================================== %
            % rigstr: a 'char' variable typically describing the...
            %   ...settings of an experiment.
            % ======================================================== %
            
            % initialize break point location for title if not input
            if nargin < 2; varargin{1} = struct(); end
            if ~isfield(varargin{1}, 'bp')
                varargin{1}.bp = 2;
            end

            % unitalicize subscripts in variables
            ii1_rstr = strfind(rigstr0, '_');
            ii2_rstr = strfind(rigstr0, '$');
            ii_rstr = zeros(numel(ii1_rstr), 2);
            rigstr = cell(numel(ii1_rstr) + 1, 1);
            for i = 1 : numel(ii1_rstr)
                if ii1_rstr(i) > ii2_rstr(2*i-1) &&...
                        ii1_rstr(i) < ii2_rstr(2*i)
                    ii_rstr(i,1) = ii1_rstr(i) - 1;
                    ii_rstr(i,2) = ii2_rstr(2*i) - 1;
                    if i == 1
                        rigstr{i} = strcat(rigstr0(1 : ii1_rstr(i)),...
                            {'\mathrm{'}, rigstr0(ii1_rstr(i) + 1 :...
                            ii2_rstr(2*i) - 1), {'}'},...
                            rigstr0(ii2_rstr(2*i)));
                    else
                        rigstr{i} = strcat(rigstr0(ii2_rstr(2*(i-1)) + 1 :...
                            ii1_rstr(i)), {'\mathrm{'},...
                            rigstr0(ii1_rstr(i) + 1 : ii2_rstr(2*i) - 1),...
                            {'}'}, rigstr0(ii2_rstr(2*i)));
                    end
                end
            end
            rigstr{end} = rigstr0(ii2_rstr(end) + 1 : end);
            rigstr = cell2mat(strcat(rigstr{:}));

            % adjust the length of description text
            i_rstr = strfind(rigstr, '$');
            n_rstr = ceil(numel(i_rstr) / varargin{1}.bp);
            rigstr = char(strcat(rigstr(1 : i_rstr(n_rstr + 1)),...
                string(newline), rigstr(i_rstr(n_rstr + 1) + 1 : end)));            
            % rigstr = char(strcat(rigstr(1 : i_rstr(1) - 1),...
            %     string(newline), rigstr(i_rstr(1) : i_rstr(n_rstr + 1)...
            %     - 1), string(newline), rigstr(i_rstr(n_rstr + 1) : end)));            
            
        end

        
        
    end

end

