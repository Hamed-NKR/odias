
% PLOTCI  Plot MAP estimate along with credible intervals.
%  
%  AUTHOR: Timothy Sipkens, 2021-04-09 - Revised by Hamed Nikookar, Dec. 24

function h = plotci_custom(d, x, Gpo, x0, cspec, linstl, opts)

% initialize plot
if ~exist('opts', 'var') 
    opts = struct();
end
if (~isfield(opts, 'log')) || isempty(opts.log)
    opts.log = 'off'; % use continuum regime approximation
end

if ~exist('x0', 'var'); x0 = []; end

% Plot color.
if ~exist('cspec', 'var'); cspec = []; end
% if isempty(cspec); cspec = [0.36, 0.79, 0.98]; end  % pastel blue
if isempty(cspec); cspec = [1, 0.37, 0.54]; end  % pastel red

% Parse error bound input.
% If not supplied, do not plot.
if ~exist('Gpo', 'var'); Gpo = []; end
if isempty(Gpo), Gpo = zeros(length(x)); end


% Credible interval limits for plotting.
x_high2 = x + 2 .* sqrt(diag(Gpo));
x_low2 = max(x - 2 .* sqrt(diag(Gpo)),0);
x_high1 = x + sqrt(diag(Gpo));
x_low1 = max(x - sqrt(diag(Gpo)),0);

if strcmp(opts.log, 'on') || strcmp(opts.log, 'On') || strcmp(opts.log, 'ON')
    ind_rmv = (x_low1 <= 0) | (x_low2 <= 0) | (x_high1 <= 0) |...
        (x_high2 <= 0);
    x_low1(ind_rmv) = [];
    x_low2(ind_rmv) = [];
    x_high1(ind_rmv) = [];
    x_high2(ind_rmv) = [];
    d(ind_rmv) = [];
    x(ind_rmv) = [];
end

% Plot shaded region.
reg = [x_low2; flipud(x_high2)];
fill([d; flipud(d)], reg, cspec, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.15);

set(gca, 'XScale', 'log');
hold on;

reg = [x_low1; flipud(x_high1)];
fill([d; flipud(d)], reg, cspec, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.22);


% Plot estimate.
h = plot(d, x, 'Color', cspec, 'LineWidth', 1.5, 'LineStyle', linstl);


% Plot credible intervals.
% plot(d, x_high2, '--', 'Color', [cm, 0.5]);
% plot(d, x_low2, '--', 'Color', [cm, 0.5]);


% Plot true solution.
if ~isempty(x0)
    plot(d, x0, 'k--', 'LineWidth', 1);
end

hold off;
axis square;

if nargout==0; clear h; end

end

