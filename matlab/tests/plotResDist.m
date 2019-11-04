% Evaluates the numerics of solveImuArray.

%% Initialization.
Na = 3;
Ng = 1;
sigmaa = 1e-6;
sigmag = 1e-6;

iters = 1000;

ress = inf(1,iters);
resw = inf(1,iters);
reswprim = inf(1,iters);

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
    y = [ya; yg];

    [s,w,wprime] = solveImuArray(ya,yg,r,sigmaa,sigmag);
    
    % Pick most likely solution.
%     s = s(:,1);
%     w = w(:,1);
%     wprime = wprime(:,1);

    ress(k) = min(vecnorm((s-sgt)));
    resw(k) = min(vecnorm((w-wgt)));
    reswprim(k) = min(vecnorm(wprime-wprimegt));

    [s,w,wprime] = solveImuArrayMl(ya,yg,r,sigmaa,sigmag,wgt);
    reswml(k) = norm(w-wgt);
end

%% Plot residual distributions.
edges = -20:0.25:5;
centers = edges(1:end-1)+diff(edges)/2;
hcs = histcounts(log10(ress),edges)/iters;
hcw = histcounts(log10(resw),edges)/iters;
hcwprime = histcounts(log10(reswprim),edges)/iters;

hcwml = histcounts(log10(reswml),edges)/iters;

figure;
% plot(centers,hcs,'-','LineWidth',2);
hold on
plot(centers,hcw,'-','LineWidth',2);
% plot(centers,hcwprime,'-','LineWidth',2);

hold on
plot(centers,hcwml,'--','LineWidth',2);
hold off
xlabel('log_{10}(residuals)');
% legend('s','\omega','\omega''');
legend('Action matrix method','Gauss-Newton from GT','Location','NorthWest');
set(gca,'FontName','Times');
set(gca,'FontSize',12);

