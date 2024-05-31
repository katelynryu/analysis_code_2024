function [x_new]=filter_position(x_raw,time,cutoff)
% Filter out outliers in the positions with a quantile threshold for
% diff(positions). The outliers are replaced by linear interpolations
%
% INPUT:
%   - x_raw: one-dimensional array of physical positions
%   - time: same size as x_raw. The time points at which each position point is recorded
%   - cutoff: float. Threshold for between-frame position displacement
%     For example, with cutoff=12.5, if the position displacement between
%   two frames are larger than 12.5 (same unit as x_raw), this
%   position will be replaced by the linear interpolation between the
%   previous non-outlier position and the next non-outlier position.
%
% OUTPUT:
%   - x_new: position vector which gets filtered once. Same size as x_raw.
%   Note that outliers may still exist.

x_new=x_raw;
x1diff=abs(diff(x_raw));
% threshold on 0.95 quantile
err_x_inx=find(x1diff>cutoff)+1;
x1_temp=x_raw;
x1_temp(err_x_inx)=[];
time_temp=time;
time_temp(err_x_inx)=[];
xq=interp1(time_temp,x1_temp,time(err_x_inx));
%update x1
x_new(err_x_inx)=xq;
end