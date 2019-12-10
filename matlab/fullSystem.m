function [s,w,wprime,R] = fullSystem(ya,yg,r,sa,sg)
% FULLSYSTEM Finds the forces for an IMU array
%   [s,w,wprime,R] = FULLSYSTEM(ya,yg,r,sa,sg) TODO: Documentation.

    Nt = size(ya,3);
    Na = size(ya,2);
    Ng = size(yg,2);
    
    nSamples = 10;

    % Find pairwise rotations between IMUs.
    % TODO: Maybe place this in the RANSAC loop and sample two time steps
    % and find the transforms that way. Should be more robust to outliers
    % but maybe less so to noise.
    Rp = cell(Ng,1);
    for i=1:Ng
        Rp{i} = wahba(squeeze(yg(:,i,:)),squeeze(yg(:,1,:)));
    end
    
    % RANSAC sample loop.
    bestSol = struct();
    bestSol.res = inf;
    for iSample=1:nSamples
        % Random minimal sample - 1 time step, 3 IMUs.
        t = randi(Nt);
        imu = randperm(Na,3);
        ygt = yg(:,imu,t);
        yat = ya(:,imu,t);
        rt = r(:,imu);
        Rpt = Rp(imu);
        
        % Solve sample.
        R1 = solveImuArrayR(yat,ygt,rt,Rpt);
        % TODO: Maybe some solutions can be discarded directly to avoid
        % solving for s, w, and wprime. Should make things faster.
        
        for iSol=1:length(R1)
            R = cell(length(Rp),1);
            for j=1:length(Rp)
                R{j} = R1{iSol}*Rp{j};
            end
            
            % Transform measurements to global frame.
            za = zeros(size(ya));
            zg = zeros(size(yg));
            for ia=1:Na
                za(:,ia,:) = R{ia}*squeeze(ya(:,ia,:));
                zg(:,ia,:) = R{ia}*squeeze(yg(:,ia,:));
            end

            % Solve s, wprime and w for all time steps.
            s = zeros(3,Nt);
            w = zeros(3,Nt);
            wprime = zeros(3,Nt);
            for it=1:Nt
                w0 = mean(zg(:,:,it),2);
                [s(:,it),w(:,it),wprime(:,it)] = solveImuArrayMl(za(:,:,it),zg(:,:,it),r,sa,sg,w0);
            end

            % One could refine the solution here but when testing it seems
            % to be a waste of time. The convergence of the refinement is
            % slow and we might just as well grab a new sample.
%             [s,w,wprime,R] = refineFullSystem(ya,yg,r,sa,sg,s,w,wprime,R,5);

            % Evaluate sample.
            res = fullSystemResidual(ya,yg,r,sa,sg,s,w,wprime,R);
            
            if res < bestSol.res
                bestSol.res = res;
                bestSol.R = R;
                bestSol.s = s;
                bestSol.w = w;
                bestSol.wprime = wprime;
            end
            fprintf('Best residual so far: %f\n',bestSol.res);
        end
    end
    
    % Refine solution.
    R = bestSol.R;
    s = bestSol.s;
    w = bestSol.w;
    wprime = bestSol.wprime;
    [s,w,wprime,R] = refineFullSystem(ya,yg,r,sa,sg,s,w,wprime,R,10);
    
    res = fullSystemResidual(ya,yg,r,sa,sg,s,w,wprime,R);
    fprintf('Final residual: %f\n',res);
end

