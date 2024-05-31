function [vel] = getAngularVelocity(vecs)
% calculate angular velocity (change of direction of vectors)
% input:
% - vecs: nframes x 2
% output:
% - vel: velocity in angle/frame (same size as tracks)

vel = zeros(size(vecs,1),1);
for i = 1:(length(vecs)-1)
    vec1 = vecs(i,:);
    vec2 = vecs(i+1,:);
    cp = cross([vec1 0],[vec2 0]);
    if cp(3) <= 0
        vel(i+1) = acos(dot(vec1,vec2)/norm(vec1)/norm(vec2)); % in [0,pi]
    else
        vel(i+1) = -acos(dot(vec1,vec2)/norm(vec1)/norm(vec2)); % in [0,pi]
    end
end
vel(1) = vel(2);
end
