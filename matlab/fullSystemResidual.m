function res = fullSystemResidual(ya,yg,r,sa,sg,s,w,wprime,R)

    Nt = size(ya,3);
    Na = size(ya,2);

    % Transform measurements to global frame.
    za = zeros(size(ya));
    zg = zeros(size(yg));
    for ia=1:Na
        za(:,ia,:) = R{ia}*squeeze(ya(:,ia,:));
        zg(:,ia,:) = R{ia}*squeeze(yg(:,ia,:));
    end
    
    % Caclulate residual.
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
end

