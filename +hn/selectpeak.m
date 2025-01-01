function dist = selectpeak(dist)

% a function to select the highest peak among the extrema of a distribution
    
    n = size(dist,1); % number of distributions existing in the dataset
    
    for i = 1 : n
        
        % store all the extrema first
        dist(i).d_mode_0 = dist(i).d_mode;
        
        % find the indices of the extrema
        d_modes = dist(i).d_mode;
        ii0 = find(ismember(dist(i).d, d_modes));
        
        % find dn/dlog(dm) corresponding to the extrema
        xs = dist(i).x(ii0);
        
        % find the peak among the extrema
        ii = find(dist(i).x == max(xs));
        dist(i).d_mode = dist(i).d(ii);

    end

end