% Plot figures 1 and 2 in manuscript

% load data
K3 = load('11thNovExpts.mat','dxvec','tocic','itic'); % these contain P_D with (1,1) block with diag(D)
K1 = load('11Novk1.mat','tocic','itic');
PD = load('11NovDtype2PD.mat','tocdiagk1','tocdiagk3','itdiagk1','itdiagk3'); % P_D with (1,1) block precond with 
% RR D with delta = 0.01
dxvec = K3.dxvec;
%% latex interpreters

set(0,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%%
figure(1);clf % 
subplot(2,2,1)

loglog(dxvec,PD.tocdiagk1(:,1),'m--v')
hold on
loglog(dxvec,PD.tocdiagk1(:,2),'r-+')
loglog(dxvec,PD.tocdiagk1(:,3),'b-s')
loglog(dxvec,PD.tocdiagk1(:,4),'k-x')
loglog(dxvec,PD.tocdiagk1(:,5),'c--o')
xlabel({'Spatial discretisation';'(a)'})
ylabel('Wallclock time')
ylim([0.15,300])
title('$$\mathbf{L}_0$$')
legend('$$\mathbf{R}_{diag}$$','$$\mathbf{R}_{block}$$','$$\mathbf{R}_{RR}$$','$$\mathbf{R}_{ME}$$','$$\mathbf{R}$$')
set(gca,'FontSize',16);

subplot(2,2,2)
% Lapp = LM wallclock
loglog(dxvec,PD.tocdiagk3(:,1),'m--v')
hold on
loglog(dxvec,PD.tocdiagk3(:,2),'r-+')
loglog(dxvec,PD.tocdiagk3(:,3),'b-s')
loglog(dxvec,PD.tocdiagk3(:,4),'k-x')
loglog(dxvec,PD.tocdiagk3(:,5),'c--o')
xlabel({'Spatial discretisation';'(b)'})
ylabel('Wallclock time')
title('$$\mathbf{L}_M$$, $$k=3$$')
ylim([0.15,300])
%legend('R_{diag}','R_{block} ichol','R_{block} chol','R_{RR} ichol','R_{RR} chol','R_{ME} ichol','R_{ME} chol','R ichol','R chol')
set(gca,'FontSize',16);

subplot(2,2,3)
% Lapp = L0 iterno
semilogx(dxvec,PD.itdiagk1(:,1),'m--v')
hold on
semilogx(dxvec,PD.itdiagk1(:,2),'r-+')
semilogx(dxvec,PD.itdiagk1(:,3),'b-s')
semilogx(dxvec,PD.itdiagk1(:,4),'k-x')
semilogx(dxvec,PD.itdiagk1(:,5),'c--o')
ylim([0,510])
xlabel({'Spatial discretisation';'(c)'})
ylabel('Iterations')
title('$$\mathbf{L}_0$$')
set(gca,'FontSize',16);

subplot(2,2,4)
% Lapp = LM iterno
semilogx(dxvec,PD.itdiagk3(:,1),'m--v')
hold on
semilogx(dxvec,PD.itdiagk3(:,2),'r-+')
semilogx(dxvec,PD.itdiagk3(:,3),'b-s')
semilogx(dxvec,PD.itdiagk3(:,4),'k-x')
semilogx(dxvec,PD.itdiagk3(:,5),'c--o')
ylim([0,510])
xlabel({'Spatial discretisation';'(d)'})
ylabel('Iterations')
title('$$\mathbf{L}_M$$, $$k=3$$')
set(gca,'FontSize',16);

%%
figure(2);clf % FIGURE ONE IN PAPER WITH BLOCK DIAGONAL
subplot(2,2,1)
% Lapp = L0 wallclock
loglog(dxvec,K1.tocic(:,1),'m--v')
hold on
loglog(dxvec,K1.tocic(:,2),'r-+')
loglog(dxvec,K1.tocic(:,3),'b-s')
loglog(dxvec,K1.tocic(:,4),'k-x')
loglog(dxvec,K1.tocic(:,5),'c--o')
xlabel({'Spatial discretisation';'(a)'})
ylabel('Wallclock time')
ylim([0.1,200])
title('$$\mathbf{L}_0$$')
legend('$$\mathbf{R}_{diag}$$','$$\mathbf{R}_{block}$$','$$\mathbf{R}_{RR}$$','$$\mathbf{R}_{ME}$$','$$\mathbf{R}$$')
set(gca,'FontSize',16);

subplot(2,2,2)
% Lapp = LM wallclock
loglog(dxvec,K3.tocic(:,1),'m--v')
hold on
loglog(dxvec,K3.tocic(:,2),'r-+')
loglog(dxvec,K3.tocic(:,3),'b-s')
loglog(dxvec,K3.tocic(:,4),'k-x')
loglog(dxvec,K3.tocic(:,5),'c--o')
xlabel({'Spatial discretisation';'(b)'})
ylim([0.1,200])
ylabel('Wallclock time')
title('$$\mathbf{L}_M$$, $$k=3$$')
%legend('R_{diag}','R_{block} ichol','R_{block} chol','R_{RR} ichol','R_{RR} chol','R_{ME} ichol','R_{ME} chol','R ichol','R chol')
set(gca,'FontSize',16);

subplot(2,2,3)
% Lapp = L0 iterno
semilogx(dxvec,K1.itic(:,1),'m--v')
hold on
semilogx(dxvec,K1.itic(:,2),'r-+')
semilogx(dxvec,K1.itic(:,3),'b-s')
semilogx(dxvec,K1.itic(:,4),'k-x')
semilogx(dxvec,K1.itic(:,5),'c--o')
ylim([0,150])
xlabel({'Spatial discretisation';'(c)'})
ylabel('Iterations')
title('$$\mathbf{L}_0$$')
set(gca,'FontSize',16);

subplot(2,2,4)
% Lapp = LM iterno
semilogx(dxvec,K3.itic(:,1),'m--v')
hold on
semilogx(dxvec,K3.itic(:,2),'r-+')
semilogx(dxvec,K3.itic(:,3),'b-s')
semilogx(dxvec,K3.itic(:,4),'k-x')
semilogx(dxvec,K3.itic(:,5),'c--o')
ylim([0,150])
xlabel({'Spatial discretisation';'(d)'})
ylabel('Iterations')
title('$$\mathbf{L}_M$$, $$k=3$$')
set(gca,'FontSize',16);
