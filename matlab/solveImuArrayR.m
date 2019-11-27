function [s,w,wprime,R,q] = solveImuArrayR(ya,yg,r,Rp)
% SOLVEIMUARRAYR Find the translation and rotation forces of an IMU array.
% TODO: Documentation. Similar to SOLVEIMUARRAY but also solves for the
% orientation of the IMUs provided their orientation relative to the first
% IMU.
%
% See also SOLVEIMUARRAY.

    Na = numel(ya)/3;
    Ng = numel(yg)/3;

    assert(Na >= 3,'There has to be at least three accelerometers.');
    assert(Na == Ng,'The number of accelerometers and gyroscopes must match.');
    assert(isequal(size(r),[3 Na]),'The size of r is inconsitent with ya. Expected a size of [3 %d].',Na);
    
    ya = reshape(ya,3,Na);
    yg = reshape(yg,3,Ng);

    Ha = [-skewSymmetric(r) repmat(eye(3),Na,1)];

    E = zeros(9,9);
    E([5 9 11 13 21 25 28 36 42 44 46 50]) = [-1 -1 1 1 1 1 -1 -1 1 1 -1 -1];
    Wa = kron(r',eye(3))*E;
%     Wg = [zeros(3*Ng,6) repmat(eye(3),Ng,1)];
%     W = [Wa; Wg]; % h(w) = W*m
    
    % Transform measurements to the frame of the first IMU.
    za = zeros(3,Na);
    for i=1:Na
        za(:,i) = Rp{i}*ya(:,i);
    end
    zg = zeros(3,Ng);
    for i=1:Ng
        zg(:,i) = Rp{i}*yg(:,i);
    end

    Va = kron(za',eye(3));
%     Vg = kron(zg',eye(3));
%     V = [Va; Vg];
    
    ASD = eye(3*Na)-Ha*((Ha'*Ha)\Ha');
    Z = ASD*[Va -Wa];
    Z = Z(1:3,1:15);
    
    % Solve polynomial system.
    qsols = solver_inertial_array_R([Z(:); mean(zg,2)]);
    
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
    
    % Find w and transforms.
    mq = zeros(9,size(q,2));
    w = zeros(3,size(q,2));
    R = cell(length(Rp),size(q,2));
    for i=1:size(q,2)
        R1 = quat2dcm(q(:,i)');
        mq(:,i) = R1(:);
        w(:,i) = R1*mean(zg,2);
        for j=1:length(Rp)
            R{j,i} = R1*Rp{j};
        end
    end
    
    wx = w(1,:);
    wy = w(2,:);
    wz = w(3,:);
    mw = [wx.^2; wx.*wy; wx.*wz; wy.^2; wy.*wz; wz.^2; wx; wy; wz];
    
    % Sort solutions based on likelihood.
%     L = -0.5*dot((V*mq-W*mw),P*(V*mq-W*mw));
%     [~,sind] = sort(L,'descend');
%     q  = q(:,sind);
%     mq = mq(:,sind);
%     w = w(:,sind);
%     R = R(:,sind);
    
    % Solve linearly for s and wprime.
    swprime = (Ha'*Ha)\Ha'*(Va*mq-Wa*mw);
    s = swprime(4:6,:);
    wprime = swprime(1:3,:);
end

