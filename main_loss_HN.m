clc
clear
close all

n_dm = 100;
dm_min = 5e-9;
dm_max = 1e-6;
dm = logspace(log10(dm_min), log10(dm_max), n_dm);

mu1 = 1.4e-7;
sigma1 = 1.6;
ntot1 = 6e7 / 1e6;

dn1_dlogdm = ntot1 * lognpdf(dm, log(mu1), log(sigma1));

% initialize size distribution figure
figure(1);
h1 = gcf;
h1.Position = [0, 0, 500, 500];
set(h1, 'color', 'white');

dist1 = plot(1e9 * dm, dn1_dlogdm);
hold on

% plot appearance configs
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'xscale', 'log')
xlim(1e9 * [dm_min, dm_max])
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 16)
ylabel('$\mathrm{d}n/\mathrm{dlog}(d_\mathrm{m})$ $[\#/\mathrm{cc}]$',...
    'interpreter', 'latex', 'FontSize', 16)

q = 0.5 / 60000; % flow rate [m3/s]
l = 1.88; % tube length [m]
do = 1; % tube O.D. [inch]
t = 0.065; % tube thickness
di = (do - 2*t) * 2.54e-2; % tube inner diameter [m]

etta = ones(1, n_dm);

for i = 1: n_dm
    etta(i) = calculate_loss(q, l, di, dm(i), 10);
end

dist2 = plot(1e9 * dm, etta .* dn1_dlogdm);

legend([dist1, dist2], {'Raw', 'Loss-corrected'}, 'interpreter', 'latex',...
    'Location', 'northeast', 'FontSize', 12)

% initialize penetration efficiency figure
figure(2);
h2 = gcf;
h2.Position = [0, 0, 500, 500];
set(h2, 'color', 'white');

plot(1e9 * dm, 100 * etta);
% plot appearance configs
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'xscale', 'log')
xlim(1e9 * [dm_min, dm_max])
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 16)
ylabel('${\eta}$ $[\%]$', 'interpreter', 'latex', 'FontSize', 16)

% calculate particle loss in a straight tube with a constant particle...
    % ...concentration
function etta = calculate_loss(q, l, d, d_p, n_s)

    %%%%%%%%% inputs/outputs %%%%%%%%%
    % etta: penetration efficiency
    % q: air flow rate
    % l: straight tube length
    % d: tube diameter
    % d_p: particle diameter
    % n_s: last element number to be considered in Gormley and Kennedy's...
    %   ...series
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    etta = 18.37e-6; % air dynamic viscosity (Pa.s)
    T = 25 + 273.15; % room temperature (K)
    k_B = 1.380649 * 10^(-23); % Boltzmann's constant (m2.kg.s-2.K-1)
    
    % Stokes–Einstein–Sutherland diffusion equation
    D = (k_B * T) / (6 * pi * etta * (d_p/2));
    
    u = q / (pi * d^2 / 4); % flow velocity

    Pe = (d * u) / D; % Peclet number
    
    x = 2 * (l/d) * Pe^(-1);

    % Gormley and Kennedy (1949)'s solution
    n = 0 : n_s;
    lambda_n = zeros(length(n),1);
    lambda_n(1) = sqrt(7.312);
    lambda_n(2) = sqrt(44.62);
    lambda_n(3) = sqrt(113.8);
    lambda_n(4:end) = 4*n(4:end) + 8/3;
    G_n = zeros(length(n),1);
    G_n(1) = 0.749;
    G_n(2) = 0.544;
    G_n(3) = 0.463;
    G_n(4:end) = 1.01276 * lambda_n(4:end).^(-1/3);

    etta = 8 * sum((G_n ./ lambda_n.^2) .* exp(-lambda_n.^2 * x));

end