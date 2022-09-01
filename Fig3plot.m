% Create figure 3: wallclock and iterno for large values of N (reduced
% Tend)
%load('AprilHighDimExample')
%%


load('12NovHighDimExp.mat');
%% latex interpreters

set(0,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure(10);clf
subplot(2,2,1)
% Lapp = L0 wallclock
dxvec = [1e-4,4e-5,2e-5];

loglog(dxvec,tocick1p4(1:3,2),'r-x')
hold on

loglog(dxvec,tocick1p4(1:3,3),'b-x')

loglog(dxvec,tocick1p4(1:3,5),'k-x')

loglog(dxvec,tocick1(1:3,2),'r--+')
loglog(dxvec,tocick1(1:3,3),'b--+')
loglog(dxvec,tocick1(1:3,5),'k--+')
ylim([3,560])
xlabel({'Spatial discretisation';'(a)'})
ylabel('Wallclock time')
title('$$\mathbf{L}_0$$')
legend('$$\mathbf{R}_{block}$$, $$s=4p$$','$$\mathbf{R}_{RR}$$, $$s=4p$$','$$\mathbf{R}$$, $$s=4p$$','$$\mathbf{R}_{block}$$, $$s=2p$$','$$\mathbf{R}_{RR}$$, $$s=2p$$','$$\mathbf{R}$$, $$s=2p$$')
set(gca,'FontSize',16);

subplot(2,2,2)
% Lapp = LM wallclock


loglog(dxvec,tocick3p4(1:3,2),'r-x')
hold on

loglog(dxvec,tocick3p4(1:3,3),'b-x')

loglog(dxvec,tocick3p4(1:3,5),'k-x')

loglog(dxvec,tocick3(1:3,2),'r--+')
loglog(dxvec,tocick3(1:3,3),'b--+')
loglog(dxvec,tocick3(1:3,5),'k--+')
ylim([3,560])
xlabel({'Spatial discretisation';'(b)'})
ylabel('Wallclock time')
title('$$\mathbf{L}_M$$')
%legend('$$\mathbf{R}_{block}$$, $$s=4p$$','$$\mathbf{R}_{RR}$$, $$s=4p$$','$$\mathbf{R}$$, $$s=4p$$','$$\mathbf{R}_{block}$$, $$s=2p$$','$$\mathbf{R}_{RR}$$, $$s=2p$$','$$\mathbf{R}$$ $$s=2p$$')
set(gca,'FontSize',16);


subplot(2,2,3)
% Lapp = L0 iterno


semilogx(dxvec,itick1p4(1:3,2),'r-x')
hold on

semilogx(dxvec,itick1p4(1:3,3),'b-x')

semilogx(dxvec,itick1p4(1:3,5),'k-x')

semilogx(dxvec,itick1(1:3,2),'r--+')
semilogx(dxvec,itick1(1:3,3),'b--+')
semilogx(dxvec,itick1(1:3,5),'k--+')
ylim([40,125])
xlabel({'Spatial discretisation';'(c)'})
ylabel('Iterations')
title('$$\mathbf{L}_0$$')
%legend('$$\mathbf{R}_{block}$$, $$s=4p$$','$$\mathbf{R}_{RR}$$, $$s=4p$$','$$\mathbf{R}$$, $$s=4p$$','$$\mathbf{R}_{block}$$, $$s=2p$$','$$\mathbf{R}_{RR}$$, $$s=2p$$','$$\mathbf{R}$$ $$s=2p$$')
set(gca,'FontSize',16);


subplot(2,2,4)


% Lapp = LM iterno

semilogx(dxvec,itick3p4(1:3,2),'r-x')
hold on

semilogx(dxvec,itick3p4(1:3,3),'b-x')

semilogx(dxvec,itick3p4(1:3,5),'k-x')

semilogx(dxvec,itick3(1:3,2),'r--+')
semilogx(dxvec,itick3(1:3,3),'b--+')
semilogx(dxvec,itick3(1:3,5),'k--+')
ylim([40,125])
xlabel({'Spatial discretisation';'(d)'})
ylabel('Iterations')
title('$$\mathbf{L}_M$$')
%legend('$$\mathbf{R}_{block}$$, $$s=4p$$','$$\mathbf{R}_{RR}$$, $$s=4p$$','$$\mathbf{R}$$, $$s=4p$$','$$\mathbf{R}_{block}$$, $$s=2p$$','$$\mathbf{R}_{RR}$$, $$s=2p$$','$$\mathbf{R}$$ $$s=2p$$')
set(gca,'FontSize',16);