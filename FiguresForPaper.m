load('12thNovCorrectedExp.mat')
set(0,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%%

figure(3);clf
subplot(1,2,1)
%plot([1,15],[itic(1,1),itic(1,1)],'m--')

plot(itic(1:16,1),'m')
hold on
plot(itic(:,2),'r')
plot(itic(:,3),'b')
plot(itic(:,4),'k')
plot(itic(:,5),'c')
legend('$$R_{diag}$$','$$R_{bl}$$','$$R_{RR}$$','$$R_{ME}$$','$$R$$')

ylabel('Iterations')
xlabel('k')
set(gca,'FontSize',16);

subplot(1,2,2)
plot(tocic(:,1),'m')
hold on
plot(tocic(:,2),'r')
plot(tocic(:,3),'b')
plot(tocic(:,4),'k')
plot(tocic(:,5),'c')
ylabel('Wallclock time')
xlabel('k')
set(gca,'FontSize',16);


figure(5);clf
subplot(1,2,1)
plot(itdiagD2(:,1)','m')
hold on
plot(itdiagD2(:,2),'r')
plot(itdiagD2(:,3),'b')
plot(itdiagD2(:,4),'k')
plot(itdiagD2(:,5),'c')
ylabel('Iterations')
xlabel('k')
legend('$$R_{diag}$$','$$R_{bl}$$','$$R_{RR}$$','$$R_{ME}$$','$$R$$')
set(gca,'FontSize',16);

subplot(1,2,2)
plot(tocdiagD2(:,1),'m')
hold on
plot(tocdiagD2(:,2),'r')
plot(tocdiagD2(:,3),'b')
plot(tocdiagD2(:,4),'k')
plot(tocdiagD2(:,5),'c')
ylabel('Wallclock time')
xlabel('k')
set(gca,'FontSize',16);
