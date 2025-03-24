function [dn_dlogd, d] = kde(d0, dn0_dlogd, lambda, rd, n_rsl)
% kde performs Kernel Density Estimation on a group size distributions...
    % ...to output smoother curves.

% d0: a cell array of size setpoints from different distributions
% dn0_dlogd: a cell array of density distributions of counts in...
    % ...log space that correspond to d0
% d: high resolution size setpoints for output data
% dn_dlogd: smoothed size distributions
% lambda: smoothing parameter
% rd: extrapolation factor
% n_rsl: resolution

% define the range for interpolation (extrapolate beyond observed range)
dmin = min(d0) / rd; % extend lower range
dmax = max(d0) * rd;   % extend upper range
log_d_construct = linspace(log10(dmin), log10(dmax), n_rsl); % defined output log-space...
    % ...increments for size for better resolution

% kernel density estimation
[pdf_smooth, log_d] = ksdensity(log10(d0), log_d_construct,...
    'Bandwidth', lambda*range(log10(d0)), 'Support', 'positive');
d = 10.^log_d;

% normalize KDE to match input PDF scaling
dn_dlogd = pdf_smooth * trapz(d0, dn0_dlogd) / trapz(d, pdf_smooth);

end

