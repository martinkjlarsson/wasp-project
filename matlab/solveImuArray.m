function [s,w,wprime] = solveImuArray(ya,yg,r,sa,sg)
% SOLVEIMUARRAY Find the translation and rotation forces of an IMU array.
%   [s,w,wprime] = solveImuArray(ya,yg,r,sa,sg) finds the common specific
%       forces s, the angular velocity w and the angular acceleration
%       wprime of an IMU array. The array must contain at least three
%       accelerometers. If there are N solutions to the problem the outputs
%       become 3 by N matrices.
%
%       Inputs:
%           ya - 3 by Na matrix of accelerometer measurements.
%           yg - 3 by Ng matrix of gyroscope measurements.
%           r - 3 by Na matrix of accelerometer positions.
%           sa - standard deviation of the accelerometer measurements.
%           sg - standard deviation of the gyroscope measurements.
%
%   See I. Skog, Inertial Sensor Arrays, Maximum Likelihood, and Cramer-Rao
%   Bound (2016) for more details.

    % Sanitate input data.
    ya = ya(:);
    yg = yg(:);
    y = [ya; yg];

    Na = length(ya)/3;
    Ng = length(yg)/3;

    assert(Na >= 3,'There has to be at least three accelerometers.');
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
    
    % Remove missing and complex solutions.
    w = zeros(3,size(wsols,2));
    keep = false(1,size(wsols,2));
    for i=1:size(wsols,2)
        wsol = wsols(:,i);
        
        if ~all(isfinite(wsol))
            continue;
        end

        if any(abs(imag(wsol)./real(wsol)) > 1e-6)
            continue;
        end
        
        w(:,i) = real(wsol);
        keep(i) = 1;
    end
    w = w(:,keep);
    
    wx = w(1,:);
    wy = w(2,:);
    wz = w(3,:);
    m = [wx.^2; wx.*wy; wx.*wz; wy.^2; wy.*wz; wz.^2; wx; wy; wz];
    
    % Sort solutions based on likelihood.
    L = -0.5*dot((y-W*m),P*(y-W*m));
    [~,sind] = sort(L,'descend');
    w = w(:,sind);
    
    % Solve linearly for s and wprime.
    swprime = HQ*(ya-Wa*m);
    s = swprime(4:6,:);
    wprime = swprime(1:3,:);
end

