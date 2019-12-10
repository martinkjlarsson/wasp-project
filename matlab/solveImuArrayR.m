function R = solveImuArrayR(ya,yg,r,Rp)
% SOLVEIMUARRAYR Find rotation from a common frame to a global
%   R = SOLVEIMUARRAYR(ya,yg,r,Rp)
%       Inputs:
%           ya - accelerometer measurements in local IMU frame
%           yg - gyroscope measurements in local IMU frame
%           r - 3 x N matrix of IMU positions in global frame
%           Rp - cell array of rotation matrices that takes the
%               measurements from their local frames to some common one.
%       Outputs:
%           R - cell array of rotation matrices that takes the measurements
%               from the common frame to the global one. Each matrix is a
%               possible solution to problem being solved.
%
% See also SOLVEIMUARRAY.

    Na = numel(ya)/3;
    Ng = numel(yg)/3;

    assert(Na >= 3,'There has to be at least three accelerometers.');
    assert(Na == Ng,'The number of accelerometers and gyroscopes must match.');
    assert(isequal(size(r),[3 Na]),'The size of r is inconsitent with ya. Expected a size of [3 %d].',Na);
    
    ya = reshape(ya,3,Na);
    yg = reshape(yg,3,Ng);
    
    % Some configurations of r might result in bad numerics. This is just a
    % quirk of the action matrix solvers. We mitigate this problem by
    % applying a random rotation to r and revert the rotation at the end.
    randRot = quat2dcm(randn(1,4));
    r = randRot*r;
    
    % Transform measurements to the frame of the first IMU.
    za = zeros(3,Na);
    for i=1:Na
        za(:,i) = Rp{i}*ya(:,i);
    end
    zg = zeros(3,Ng);
    for i=1:Ng
        zg(:,i) = Rp{i}*yg(:,i);
    end
    mzg = mean(zg,2);

    E = zeros(9,6);
    E([5 9 11 13 21 25 28 36 42 44 46 50]) = [-1 -1 1 1 1 1 -1 -1 1 1 -1 -1];
    Wa = kron(r',eye(3))*E;
    Va = kron(za',eye(3));

    Ha = [-skewSymmetric(r) repmat(eye(3),Na,1)];
    P = eye(3*Na)-Ha*((Ha'*Ha)\Ha');
    Z = P*[Va -Wa];
    
    % Why rows 1, 5 and 9? When Na == 3 the problem is not well-defined for
    % any three equations. These happen to work. I have not investigated
    % further.
    Z = Z([1 5 9],:);

    % Solve polynomial system. Although the equations are the same there
    % are 32 solutions when Na == 3 and 48 solutions when Na >= 4. Because
    % of this, two solvers has been created and they do not perform good
    % numerically on each other's data, i.e., both are needed.
    % Consider using MEX versions for significant speed improvement.
    if Na == 3
        qsols = solver_inertial_array_R_3([Z(:); mzg]);
%         qsols = solver_inertial_array_R_3_mex([Z(:); mzg]);
    else
        qsols = solver_inertial_array_R([Z(:); mzg]);
%         qsols = solver_inertial_array_R_mex([Z(:); mzg]);
    end
    
    % Remove missing and complex solutions.
    q = zeros(4,size(qsols,2));
    keep = false(1,size(qsols,2));
    for i=1:size(qsols,2)
        qsol = qsols(:,i);
        
        if ~all(isfinite(qsol))
            continue;
        end

        if any(abs(imag(qsol)./real(qsol)) > 1e-6)
            continue;
        end
        
        q(:,i) = real(qsol);
        keep(i) = 1;
    end
    q = q(:,keep);
    q = q./vecnorm(q);
    
    % Calculate final rotation matrices by reverting the random rotation.
    R = cell(1,size(q,2));
    for i=1:size(q,2)
        R{i} = randRot'*quat2dcm(q(:,i).');
    end
end

