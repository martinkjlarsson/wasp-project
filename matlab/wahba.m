function R = wahba(p1,p2,w)
% WAHBA Finds rotation between two sets of points
%   R = WAHBA(p1,p2) finds the rotation matrix R such that p2 = R*p1.
%   R = WAHBA(p1,p2,w) applies a weighting to each point correspondence.

    if nargin < 3
        w = ones(1,size(p1,2));
    end

    W = diag(w);

    [U,~,V] = svd(p1*W*p2');
    S = eye(3);
    if det(V*U') < 0
        S(3,3) = -1;
    end
    R = V*S*U';
end

