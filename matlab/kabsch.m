function R = kabsch(p1,p2)
% KABSCH Finds rotation between two sets of points
%   R = KABSCH(p1,p2) finds the rotation matrix R such that p2 = R*p1.

    [U,~,V] = svd(p1*p2');
    S = eye(3);
    if det(V*U') < 0
        S(3,3) = -1;
    end
    R = V*S*U';
end

