function [s,w,wprime] = solveImuArray(ya,yg,r,sa,sg)
% SOLVEIMUARRAY Find the translation and rotation forces of an IMU array.
%   [s,w,wprime] = solveImuArray(ya,yg,r,sa,sg) finds the common specific
%       forces s, the angular velocity w and angular acceleration wprim of
%       and IMU array. The array must contain at least three accelerometers
%       and one gyroscope.
%
%       Inputs:
%           ya - 3 by Na matrix of accelerometer measurements.
%           yg - 3 by Ng matrix of gyroscope measurements.
%           r - 3 by Na matrix of accelerometer positions.
%           sa - standard deviation of the accelerometer measurements.
%           sg - standard deviation of the gyroscope measurements.
%
%       See I. Skog, Inertial Sensor Arrays, Maximum Likelihood, and
%       Cramer-Rao Bound (2016) for more details.

    % Sanitate input data.
    ya = ya(:);
    yg = yg(:);
    y = [ya; yg];

    Na = length(ya)/3;
    Ng = length(yg)/3;

    assert(Na >= 3,'There has to be at least three accelerometers.');
    assert(Ng > 0,'There has to be at least one gyroscope.');
    assert(isequal(size(r),[3 Na]),'The size of r is inconsitent with ya. Expected a size of [3 %d].',Na);

    % Create problem matrices.
    Ha = [-skewSymmetric(r) repmat(eye(3),Na,1)];

    Qai = (1/sa)*eye(3*Na);
    Qgi = (1/sg)*eye(3*Ng);

    HQ = (Ha'*Qai*Ha)\Ha'*Qai;
    Pa = Qai-Qai*Ha*HQ;
    Pg = Qgi;
    
    P = blkdiag(Pa,Pg);

    E = zeros(9,9);
    E([5 9 11 13 21 25 28 36 42 44 46 50]) = [-1 -1 1 1 1 1 -1 -1 1 1 -1 -1];
    Wa = kron(r',eye(3))*E;
    Wg = [zeros(3*Ng,6) repmat(eye(3),Ng,1)];
    W = [Wa; Wg]; % h(w) = W*m
    
    Z1 = W'*P*W;
    Z2 = W'*P*y;

    % Note that Z1 and W'*P could be precalculated as long as the
    % measurement data y = [ya; yg] is the only thing that changes.
    
    % Solve polynomial system.
    indsZ1 = [1;10;11;19;20;21;28;29;30;31;37;38;39;40;41;46;47;48;49;50;51;81];
    wsols = solver_inertial_array([Z1(indsZ1); Z2]);
    
    bestResNorm = inf;
    bestInd = 0;
    
    % Find best solution.
    for i=1:size(wsols,2)
        if ~isreal(wsols(:,i))
            continue;
        end
        
        if all(wsols(:,i)==0)
            continue;
        end

        wx = wsols(1,i);
        wy = wsols(2,i);
        wz = wsols(3,i);

        m = [wx^2; wx*wy; wx*wz; wy^2; wy*wz; wz^2; wx; wy; wz];
        mprime = [2*wx 0 0; wy wx 0; wz 0 wx; 0 2*wy 0; 0 wz wy; 0 0 2*wz; 1 0 0; 0 1 0; 0 0 1];
        resNorm = norm(mprime'*(Z2-Z1*m));
        if resNorm < bestResNorm
            bestResNorm = resNorm;
            bestInd = i;
        end
    end
    
    % Solve linearly for s and wprim.
    if bestInd == 0 
        s = nan;
        w = nan;
        wprime = nan;
    else
        w = wsols(:,bestInd);
        wx = w(1);
        wy = w(2);
        wz = w(3);
        m = [wx^2; wx*wy; wx*wz; wy^2; wy*wz; wz^2; wx; wy; wz];
        swprime = HQ*(ya-Wa*m);
        s = swprime(4:6);
        wprime = swprime(1:3);
    end
end

