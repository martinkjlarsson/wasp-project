% Evaluates the numerics of solveImuArray over different noise levels.

%% Initialization.
Na = 3;
Ng = 0;
nSigmas = 10;
sigmaas = linspace(1e-6,1e-1,nSigmas);
sigmags = linspace(1e-6,1e-1,nSigmas);

iters = 1000;

resw = inf(nSigmas,iters);
reswml = inf(nSigmas,iters);
sols = [];
%% Eval numerics.
for i=1:nSigmas
    fprintf('Working on sigma %d/%d\n',i,nSigmas);
    
    sigmaa = sigmaas(i);
    sigmag = sigmags(i);
    
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
        resw(i,k) = min(vecnorm((w-wgt)));

        [~,w,~] = solveImuArrayMl(ya,yg,r,sigmaa,sigmag,wgt);
        if Ng == 0
            w = [w -w]; % Negative angular velocity is also a possible solution.
        end
        reswml(i,k) = min(vecnorm((w-wgt)));
    end
end

%% Plot residuals over noise.
meanResw = median(resw,2);
meanReswml = median(reswml,2);

plot(sigmaas,meanResw,'-','LineWidth',2);
hold on
plot(sigmaas,meanReswml,'--','LineWidth',2);
hold off

xlabel('\sigma_a = \sigma_g');
ylabel('Error in \omega');
legend('Action matrix method','Gauss-Newton from GT','Location','NorthWest');
set(gca,'FontName','Times');
set(gca,'FontSize',12);


