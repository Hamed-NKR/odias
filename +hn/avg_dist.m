function [dn_dlogd, d, se_dn_dlogd] = avg_dist(d0, dn0_dlogd, opts)
% avg_dist takes average of a group size distributions.

% d0: a cell array of size setpoints from different distributions
% dn0_dlogd: a cell array of density distributions of counts in...
    % ...log space that correspond to d0
% d: common size setpoints for average output distribution
% dn_dlogd: average size distribution
% se_dn_dlogd: standard errors for average size distribution

% initialize the options structure
if ~exist('opts', 'var'); opts = struct(); end

% set default method to interpolate data using high order curve fitting
if (~isfield(opts, 'mtd')) || isempty(opts.mtd)
    opts.mtd = 'spline';
end

% set default resolution for the size distribution
if (~isfield(opts, 'rsl')) || isempty(opts.rsl)
    opts.rsl = 500;
end

% define a common vector of sizes as baseline for the output average...
    % ...distribution
d_min = min(cellfun(@min, d0));
d_max = max(cellfun(@max, d0));
d = logspace(log10(d_min), log10(d_max), opts.rsl);

% interpolate distribution to the common size vector
dn00_dlogd = zeros(length(d0), numel(d));
for i = 1 : length(d0)
    if strcmp(opts.mtd, 'spline')
        dn00_dlogd(i, :) = interp1(log10(d0{i}), dn0_dlogd{i}, log10(d),...
            'spline', 'extrap');
    else
        dn00_dlogd(i, :) = interp1(log10(d0{i}), dn0_dlogd{i}, log10(d),...
            'linear', 'extrap');
    end
end

% compute the mean and standard error over all size distributions 
dn_dlogd = mean(dn00_dlogd, 1);
se_dn_dlogd = std(dn00_dlogd, 0, 1) / sqrt(length(d0));

end

