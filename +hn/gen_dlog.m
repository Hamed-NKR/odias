function [dlog_dp, dp0] = gen_dlog(dp)
% A function to generate dlog(dp) from dp

n = length(dp);

dp0 = zeros(n-1,1);
dlog_dp = zeros(n,1);

for i = 1 : n-1
    dp0(i) = sqrt(dp(i) * dp(i+1));
end

dp0 = [dp(1)^2/dp0(1); dp0; dp(end)^2/dp0(end)];

for i = 1 : n
    dlog_dp(i) = log(dp0(i+1)/dp0(i));
end

end
