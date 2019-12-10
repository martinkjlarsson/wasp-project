% Evaluates the numerics of fullSystem.

%% Initialization.
Nt = 10;
Na = 3;
Ng = Na;
sigmaa = 1e-6;
sigmag = 1e-6;

%% Create data.
r = rand(3,Na);
% Some difficult configurations of r. Will no longer cause bad numerics.
% r = [zeros(3,1) eye(3)];
% r = eye(3);
% r = [0 0 0; 0 1 0; 0 0 1];
% theta = [-pi/3 pi/2 -2*pi/3];
% r = [cos(theta); sin(theta); zeros(1,3)];

% Generate ground truth unkowns.
sgt = rand(3,Nt);
wgt = rand(3,Nt);
wprimegt = rand(3,Nt);

qgt = randn(4,Na);
qgt = qgt./vecnorm(qgt);

Rgt = cell(Na,1);
for i=1:Na
    Rgt{i} = quat2dcm(qgt(:,i)');
end

% Generate measurements.
ya = zeros(3,Na,Nt);
yg = zeros(3,Ng,Nt);
for it=1:Nt
    Ow = skewSymmetric(wgt(:,it));
    Owprime = skewSymmetric(wprimegt(:,it));
    W = Ow*Ow+Owprime;

    ya(:,:,it) = sgt(:,it)+W*r+sigmaa*randn(3,Na);
    yg(:,:,it) = repmat(wgt(:,it),1,Ng)+sigmag*randn(3,Ng);
end

% Transform measurements to local frames.
for i=1:Na
    ya(:,i,:) = Rgt{i}'*squeeze(ya(:,i,:));
    yg(:,i,:) = Rgt{i}'*squeeze(yg(:,i,:));
end

%% Solve system.
[s,w,wprime,R] = fullSystem(ya,yg,r,sigmaa,sigmag);
q = zeros(4,length(R));
for i=1:length(R)
    q(:,i) = dcm2quat(R{i});
end

%% Print errors
% Remove quaternion sign ambiguity.
q = q.*sign(q(1,:));
qgt = qgt.*sign(qgt(1,:));

ress = rms(sgt(:)-s(:));
resw = rms(wgt(:)-w(:));
reswprime = rms(wprimegt(:)-wprime(:));
resq = rms(qgt(:)-q(:));

fprintf('\nRMS errors\ns:\t\t%e\nw:\t\t%e\nwprime:\t%e\nq:\t\t%e\n',...
    ress,resw,reswprime,resq);

