function dist = import_dist(fadd, sid, df, dm_rowid)
% A function to import tsi-generated dn/dlog(dp) distributions

% fadd: address of tsi file (including folder and file names)
% sid: smps scans to be read
% df: dilution factor
% dm_rowid: the line number in the file that has mobility setpoints
% dist: a structure containing distribution data read from tsi file

% import TSI inverted data
opts = detectImportOptions(fadd);
opts = setvartype(opts, 'double');
dist.tab0 = readtable(fadd, opts);

% save file address
dist.fadd = fadd;

% locate mobility setpoint data if not given by the user 
if ~exist('dm_rowid', 'var') || isempty(dm_rowid)
    dm_rowid = 19;
end

% import mobility setpoints and inverted dn/dlog(dm)
dist.dm = table2array(dist.tab0(dm_rowid, 39:end));
dist.dn0_dlogdm = table2array(dist.tab0((dm_rowid+1):end,...
    39:end));

% extract selected scans and apply dilution factor
if ~exist('sid', 'var') || isempty(sid)
    sid = 1 : size(dist.dn0_dlogdm,1);
end
if ~exist('df', 'var') || isempty(df)
    df = 1;
end
dist.dn0_dlogdm = df * dist.dn0_dlogdm(sid,:);

% get average and SD of dn/dlog(dm)
dist.dn_dlogdm = mean(dist.dn0_dlogdm);
dist.sigma = max(std(dist.dn0_dlogdm),...
    1e-3 * max(std(dist.dn0_dlogdm)));

% get total counts
dist.n0_tot = table2array(dist.tab0((dm_rowid+1):end, 36));
dist.n0_tot = df * dist.n0_tot(sid,:);
dist.n_tot = mean(dist.n0_tot);

% area below the size distribution curve
dist.A_tot = trapz(log10(dist.dm), dist.dn_dlogdm);

end
