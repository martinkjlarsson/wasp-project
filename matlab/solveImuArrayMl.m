function [s,w,wprime] = solveImuArrayMl(ya,yg,r,sa,sg,w0)
% SOLVEIMUARRAYML Find the translation and rotation forces of an IMU array.
%   [s,w,wprime] = solveImuArray(ya,yg,r,sa,sg) finds the common specific
%       forces s, the angular velocity w and the angular acceleration
%       wprime of and IMU array. The array must contain at least three
%       accelerometers and one gyroscope.
%
%       Inputs:
%           ya - 3 by Na matrix of accelerometer measurements.
%           yg - 3 by Ng matrix of gyroscope measurements.
%           r - 3 by Na matrix of accelerometer positions.
%           sa - standard deviation of the accelerometer measurements.
%           sg - standard deviation of the gyroscope measurements.
%           w0 - starting point for optimization.
%
%   See also solveImuArray.

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
    W = [Wa; Wg]; % h(w) = W*m)

    w = w0;
    maxIters = 10;
    
    for i=1:maxIters
        wx = w(1);
        wy = w(2);
        wz = w(3);
        m = [wx^2; wx*wy; wx*wz; wy^2; wy*wz; wz^2; wx; wy; wz];
        mprime = [2*wx 0 0; wy wx 0; wz 0 wx; 0 2*wy 0; 0 wz wy; 0 0 2*wz; 1 0 0; 0 1 0; 0 0 1];
        Jh = W*mprime;
        w = w+(Jh'*P*Jh)\(Jh'*P*(y-W*m));
        
        % TODO: Add exit condition.
        
%         L = -0.5*(y-W*m)'*P*(y-W*m);
%         fprintf('Iteration: %d, residual norm: %e, log-likelihood: %e\n',i,norm(y-W*m),L);
    end
    
    wx = w(1);
    wy = w(2);
    wz = w(3);
    m = [wx^2; wx*wy; wx*wz; wy^2; wy*wz; wz^2; wx; wy; wz];
    swprime = HQ*(ya-Wa*m);
    s = swprime(4:6);
    wprime = swprime(1:3);
end

