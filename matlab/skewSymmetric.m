function X = skewSymmetric(x)
% SKEWSYMMETRIC Gets the skew-symmetric cross product matrix
%   X = skewSymmetric(x) gets the skew-symmetric matrix corresponding to
%       cross product with x, i.e., the X such that X*a = cross(x,a). If x
%       is a 3 by n matrix X will be a 3n by 3 matrix of column stacked
%       skew-symmetric matrices.

    X = x(1)*zeros(3*size(x,2),3);
    for i=1:size(x,2)
        X(3*i-2:3*i,:) =...
            [0 -x(3,i) x(2,i); x(3,i) 0 -x(1,i); -x(2,i) x(1,i) 0];
    end
end