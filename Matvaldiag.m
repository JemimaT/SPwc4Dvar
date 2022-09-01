% number of evaluations for R, B/Q, Rinv, B/Qinv, M in forward and M in
% precond

%[R,Rinv,D,Dinv,M,Mpre] = Matvaldiagout(100,6,4);
%[R2,Rinv2,D2,Dinv2,M2,Mpre2] = Matvaldiagout(100,5,1);
clear

PD = load('11NovDtype2PD.mat','tocdiagk1','tocdiagk3','itdiagk1','itdiagk3');
K3 = load('11thNovExpts.mat','dxvec','tocic','itic'); % these contain P_D with (1,1) block with diag(D)
K1 = load('11Novk1.mat','tocic','itic');
N = 5;
% k =1;
[R,Rinv,D,Dinv,M,Mpre] = Matvaldiagout(PD.itdiagk1(:,1),5,1,0); % bad R precond

VarNames = {'R', 'Rinv', 'D', 'Dinv','M'};
T1 = table(R,Rinv,D,Dinv,M+ Mpre, 'VariableNames',VarNames)
[R,Rinv,D,Dinv,M,Mpre] = Matvaldiagout(PD.itdiagk1(:,2),5,1,1);
T1r = table(R,Rinv,D,Dinv,M+ Mpre, 'VariableNames',VarNames)

[R,Rinv,D,Dinv,M,Mpre] = Matvaldiagout(PD.itdiagk3(:,1),5,3,0); % bad R precond

VarNames = {'R', 'Rinv', 'D', 'Dinv','M'};
T3 = table(R,Rinv,D,Dinv,M+ Mpre, 'VariableNames',VarNames)
[R,Rinv,D,Dinv,M,Mpre] = Matvaldiagout(PD.itdiagk3(:,2),5,3,1);
T3r = table(R,Rinv,D,Dinv,M+Mpre, 'VariableNames',VarNames)

%%%%%%%% IC
VarNames2 = {'R', 'Rinv', 'D', 'M'};
[R,Rinv,D,M,Mpre] = MatvalICout(K1.itic(:,1),5,1,0);
T1ic = table(R,Rinv,D,M+ Mpre, 'VariableNames',VarNames2)

[R,Rinv,D,M,Mpre] = MatvalICout(K1.itic(:,2),5,1,1);
T1icr = table(R,Rinv,D,M+ Mpre, 'VariableNames',VarNames2)

[R,Rinv,D,M,Mpre] = MatvalICout(K3.itic(:,1),5,3,0);
T3ic = table(R,Rinv,D,M+ Mpre, 'VariableNames',VarNames2)

[R,Rinv,D,M,Mpre] = MatvalICout(K3.itic(:,2),5,3,1);
T3icr = table(R,Rinv,D,M+ Mpre, 'VariableNames',VarNames2)

function [R,Rinv,D,Dinv,M,Mpre] = Matvaldiagout(iter,N,k,rtype)
D = 2*iter*(N+1);
R = iter*(N+1);
M = 2*iter*N;

% preconditioner
if rtype == 0
    Rinv = 0*iter;
else
    Rinv = iter*(N+1);
end
Dinv = iter*(N+1);

% how many blocks
nblocks = floor((N+1)/k);
rem = N+1-k*nblocks;
if rem >0
    Mpre = iter*(2*(k-1)*nblocks+2*(rem-1));
else
    Mpre = iter*(2*(k-1)*nblocks);
end
end

function [R,Rinv,D,M,Mpre] = MatvalICout(iter,N,k,rtype)
D = 2*iter*(N+1);
R = iter*(N+1);
M = 2*iter*N;

% preconditioner
if rtype == 0
    Rinv = 0*iter;
else
    Rinv = iter*(N+1);
end

% how many blocks
nblocks = floor((N+1)/k);
rem = N+1-k*nblocks;
if rem >0
    Mpre = iter*(2*(k-1)*nblocks+2*(rem-1));
else
    Mpre = iter*(2*(k-1)*nblocks);
end
end