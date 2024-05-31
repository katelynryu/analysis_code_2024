function x=repeat_filter_pos(x_raw,time,threshold_quantile,iter)
% Filter out outliers in the positions with a quantile threshold for
% diff(positions). The outliers are replaced by linear interpolations
%
% INPUT:
%   - x_raw: one-dimensional array of physical positions
%   - time: same size as x_raw. The time points at which each position point is recorded
%   - threshold_quantile: (0-1 float). Quantile threshold for between-frame position displacement
%     For example, with threshold_quantile=0.98, if the position displacement between
%   two frames are larger than 98% of the between-frame displacement, this
%   position will be replaced by the linear interpolation between the
%   previous non-outlier position and the next non-outlier position.
%   - iter: int. Max number of iterations for repeating the above filtering
%
% OUTPUT:
%   - x_new: Same size as x_raw. position vector which gets filtered until 
%           1) all diff(x) are below the original threshold from diff(x_raw)
%           2) max number of iterations are reached. Note that outliers may still exist.
% 
cutoff=quantile(abs(diff(x_raw)),threshold_quantile);
x1=x_raw;
c=0;
while (any (abs(diff(x1)) > cutoff)) && (c <= iter)
    x1=filter_position(x1,time,cutoff);
    c=c+1;
end
x=x1;
end