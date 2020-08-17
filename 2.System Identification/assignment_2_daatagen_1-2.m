clear all
close all
clc
load input.mat;
load output.mat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  1) Model definition with estimation and prediction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The ARX is defined as 
% y(t) = -a1*y(t-1)-a2*y(t-2) + b0*u(t) + e(t)

N = 4000; % length of data vector
noise_std = 0.01;
uest = u(1:N/2);
yest = y(1:N/2);

uval = u(N/2+1:end);
yval = y(N/2+1:end);

PHI= zeros(3,N/2);

%% For 1 step predictor of yhat(t) = [u(t),y(t-1),y(t-2)] [b0; -a1;-a2]
estPHI(:,1)=[uest(1); 0; 0];
estPHI(:,2)=[uest(1); yest(1); 0];
PHIv(:,1)=[uval(1); 0; 0];
PHIv(:,2)=[uval(2); yval(1); 0];
y = zeros(N,1);
for t = 3:N/2 % zero initial condition
    PHIe(:,t)=[uest(t); yest(t-1); yest(t-2)];
    PHIv(:,t)=[uval(t); yval(t-1); yval(t-2)];
end
PHIe = PHIe.';
th=(PHIe.'*PHIe)\PHIe.'*yest;
bh = th(1)
a1h = th(2)
a2h = th(3)

ypred = PHIv.'*(th);
ep = yval - ypred;
erms = rms(ep)

%%
%simulation
ysim = zeros(N/2,1);
ysim(1)= bh*uval(1);
ysim(2)= bh*uval(2) + a1h*uval(1) +  a2h*ysim(1);
for t = 3:N/2
    ysim(t)= bh*uval(t) + a1h*ysim(t-1) + a2h*ysim(t-2);
end

%simulation error
es = yval - ysim;
esrms = rms(es)

%covarience
sigma = 0.01;
cv = sigma*inv(PHIe.'*PHIe);
v=norm(cv)