% Evaluates the numerics of solveImuArray.

%% Initialization.
Na = 3;
Ng = 0;
sigmaa = 1e-3;
sigmag = 1e-3;

iters = 1000;

resw = inf(1,iters);
reswml = inf(1,iters);

%% Eval numerics.
for k=1:iters
    r = rand(3,Na);
    
    % Generate ground truth unkowns.
    sgt = rand(3,1);
    wgt = rand(3,1);
    wprimegt = rand(3,1);

    Av = skewSymmetric(wgt);
    Aa = skewSymmetric(wprimegt);

    si = sgt+Av*Av*r+Aa*r+sigmaa*randn(3,Na);

    % Create measurement vector.
    ya = si(:);
    yg = repmat(wgt,Ng,1)+sigmag*randn(3*Ng,1);

    [~,w,~] = solveImuArray(ya,yg,r,sigmaa,sigmag);
    
    % Pick most likely solution.
    s = s(:,1);
    w = w(:,1);
    wprime = wprime(:,1);

    resw(k) = min(vecnorm((w-wgt)));

    [~,w,~] = solveImuArrayMl(ya,yg,r,sigmaa,sigmag,randn(3,1));
    reswml(k) = min(vecnorm((w-wgt)));
end

%% Plot residual distributions.
edges = -20:0.25:5;
centers = edges(1:end-1)+diff(edges)/2;
hcw = histcounts(log10(resw),edges)/iters;
hcwml = histcounts(log10(reswml),edges)/iters;

figure;
plot(centers,hcw,'-','LineWidth',2);
hold on
plot(centers,hcwml,'--','LineWidth',2);
hold off

xlabel('log_{10}(residuals)');
legend('Action matrix method','Gauss-Newton from GT','Location','NorthWest');
set(gca,'FontName','Times');
set(gca,'FontSize',12);

