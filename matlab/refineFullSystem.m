function [s,w,wprime,R] = refineFullSystem(ya,yg,r,sa,sg,s,w,wprime,R,maxIters)
% REFINEFULLSYSTEM Refine a solution from fullSystem
%
% See also fullSystem.

    Nt = size(ya,3);
    Na = size(ya,2);
    Ng = size(yg,2);
    
    if nargin < 10
        maxIters = 10;
    end
    
    for i=1:maxIters
        % Calculate specific forces in the global frame.
        si = zeros(3,Na,Nt);
        for it=1:Nt
            Ow = skewSymmetric(w(:,it));
            Owprime = skewSymmetric(wprime(:,it));
            W = Ow*Ow+Owprime;
            si(:,:,it) = s(:,it)+W*r;
        end

        % Find optimal rotations.
        for ia=1:Na
            p1 = [squeeze(ya(:,ia,:)) squeeze(yg(:,ia,:))];
            p2 = [squeeze(si(:,ia,:)) w];
            weighting = [(1/sa^2)*ones(Nt,1); (1/sg^2)*ones(Nt,1)];
            R{ia} = wahba(p1,p2,weighting);
        end

        % Transform measurements to global frame.
        za = zeros(size(ya));
        for ia=1:Na
            za(:,ia,:) = R{ia}*squeeze(ya(:,ia,:));
        end
        zg = zeros(size(yg));
        for ig=1:Ng
            zg(:,ig,:) = R{ig}*squeeze(yg(:,ig,:));
        end
        
        % Evaluate sample.
        resa = zeros(1,Nt);
        resg = zeros(1,Nt);
        for it=1:Nt
            Ow = skewSymmetric(w(:,it));
            Owprime = skewSymmetric(wprime(:,it));
            W = Ow*Ow+Owprime;
            resa(it) = sum((s(:,it)+W*r-za(:,:,it)).^2,'all');
            resg(it) = sum((w(:,it)-zg(:,:,it)).^2,'all');
        end
        res = sum(resa)/sa^2+sum(resg)/sg^2;
        fprintf('Iteration %d (rotation): \tres=%f\n',i,res);

        % Solve s, wprime and w for all time steps.
        w0 = w;
        s = zeros(3,Nt);
        w = zeros(3,Nt);
        wprime = zeros(3,Nt);
        for it=1:Nt
            [s(:,it),w(:,it),wprime(:,it)] = solveImuArrayMl(za(:,:,it),zg(:,:,it),r,sa,sg,w0(:,it));
%             [s(:,it),w(:,it),wprime(:,it)] = solveImuArray(za(:,:,it),zg(:,:,it),r,sa,sg);
        end
        
        % Evaluate sample.
        resa = zeros(1,Nt);
        resg = zeros(1,Nt);
        for it=1:Nt
            Ow = skewSymmetric(w(:,it));
            Owprime = skewSymmetric(wprime(:,it));
            W = Ow*Ow+Owprime;
            resa(it) = sum((s(:,it)+W*r-za(:,:,it)).^2,'all');
            resg(it) = sum((w(:,it)-zg(:,:,it)).^2,'all');
        end
        res = sum(resa)/sa^2+sum(resg)/sg^2;
        fprintf('Iteration %d (forces): \t\tres=%f\n',i,res);
    end
end

