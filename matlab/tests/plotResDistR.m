% Evaluates the numerics of solveImuArrayR.

%% Initialization.
Na = 4;
Ng = Na;
% sigmaa = 0;
% sigmag = 0;
sigmaa = 1e-6;
sigmag = 1e-6;

iters = 1000;

resq = inf(1,iters);

%% Eval numerics.
for k=1:iters
    r = rand(3,Na);

    % Generate ground truth unkowns.
    sgt = rand(3,1);
    wgt = rand(3,1);
    wprimegt = rand(3,1);
    
    qgt = randn(4,Na);
    qgt = qgt./vecnorm(qgt);
    
    Rgt = cell(Na,1);
    for i=1:Na
        Rgt{i} = quat2dcm(qgt(:,i)');
    end

    Rp = cell(Na,1);
    for i=1:Na
        Rp{i} = Rgt{1}'*Rgt{i};
    end

    Av = skewSymmetric(wgt);
    Aa = skewSymmetric(wprimegt);

    si = sgt+Av*Av*r+Aa*r+sigmaa*randn(3,Na);

    % Create measurement vector.
    bigR = blkdiag(Rgt{:});
    ya = bigR'*si(:);
    yg = bigR'*(repmat(wgt,Ng,1)+sigmag*randn(3*Ng,1));

    [~,w,~,R,q] = solveImuArrayR(ya,yg,r,Rp);

    if ~isempty(q)
        resq(k) = min(vecnorm(([q -q]-qgt(:,1))));
    end
end

%% Plot residual distributions.
edges = -20:0.25:5;
centers = edges(1:end-1)+diff(edges)/2;
hcw = histcounts(log10(resq),edges)/iters;

figure;
plot(centers,hcw,'-','LineWidth',2);

xlabel('log_{10}(residuals)');
set(gca,'FontName','Times');
set(gca,'FontSize',12);

