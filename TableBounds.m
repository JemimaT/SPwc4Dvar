%  print information for table 2
clear
load('N250_10thNov.mat');

% Compute L bounds - exactly or using Section 2 theory
% terms in bound
N = length(Amat.B);
%B = n250.B; Q = n250.Q; %Table 6 with value of B from expt in Sec 7
B = speye(N);Q = speye(N); % Replace D with identity
% B
kB2 = max(eigs(B,1),eigs(Q,1))/min(eigs(B,1,'sr'),eigs(Q,1,'sr'))
eigs(Amat.R,1)
eigs(Amat.R,1,'sr')

% \tilde{R}^{-1}R lA, LA
R = Amat.R;
pvec = Amat.pvec;
p = length(R);
Rdiag = diag(diag(R));
Rblock = blockfun(R,pvec,1,p,length(p));
RRR = R + speye(p);
[v1,e] = eigs(R,2,'sr');
Gam = e(2,2)-e(1,1);
v = sqrt(Gam)*v1(:,1);
RME = R + v*v';

lA_diag = min([1,eigs(Rdiag\R,1,'sr')]);
lA_Rblock = min([1,eigs(Rblock\R,1,'sr')]);
lA_RR = min([1,eigs(RRR\R,1,'sr')]);
lA_RME = min([1,eigs(RME\R,1,'sr')]);
lA_R = 1;

LA_diag = max([1,eigs(Rdiag\R,1)]);
LA_Rblock = max([1,eigs(Rblock\R,1)]);
LA_RR = max([1,eigs(RRR\R,1)]);
LA_RME = max([1,eigs(RME\R,1)]);
LA_R = 1;

lA = [lA_diag, lA_Rblock, lA_RR, lA_RME,lA_R]
LA = [LA_diag, LA_Rblock, LA_RR, LA_RME,LA_R]

% \tilde{S}^{-1}S where using exact L and only dropping H^T R^{-1} H
% NB no dependence on \tilde{L} or \tilde{R}
k = Amat.nsubw+1; M= Amat.M;
C = full(gallery('tridiag',k,-1,0,0));
L = kron(speye(k),speye(N)) + kron(C,M);

LL_0 = L'*L;
Sest = LL_0\(L'*L+kron(speye(k),Amat.H'*(R\Amat.H))); % Lambda S terms for all cases
eigS = real(eig(Sest));
lS = min(eigS);
LS = max(eigS);
%%
ll0 = min(eig(LL_0))
LL0 = max(eig(LL_0)) % lambda L for k=0

%% construct L3

L3 = speye(6*N);
L3(251:500,1:250) = M;
L3(501:750,1:250) = M^2;
L3(501:750,251:500) = M;

L3(1001:1250,751:1000) = M;
L3(1251:1500,751:1000) = M^2;
L3(1251:1500,1001:1250) = M;

Lhat_pre = L3'*(L'*L)*L3;

ll3 = min(eig(Lhat_pre))
LL3 = max(eig(Lhat_pre))

Lhat_pre = (L3'*L3)


%% Compute bounds
% L0
% Rdiag
[0.5*(lA(1)-sqrt(lA(1)^2+4*LA(1)*LL0*LS)), 0.5*(LA(1)-sqrt(LA(1)^2+4*lA(1)*lS*ll0)),lA(1),0.5*(LA(1)+sqrt(LA(1)^2+4*LA(1)*LS*LL0))]
% Rblock
[0.5*(lA(2)-sqrt(lA(2)^2+4*LA(2)*LL0*LS)), 0.5*(LA(2)-sqrt(LA(2)^2+4*lA(2)*lS*ll0)),lA(2),0.5*(LA(2)+sqrt(LA(2)^2+4*LA(2)*LS*LL0))]
% RR
[0.5*(lA(3)-sqrt(lA(3)^2+4*LA(3)*LL0*LS)), 0.5*(LA(3)-sqrt(LA(3)^2+4*lA(3)*lS*ll0)),lA(3),0.5*(LA(3)+sqrt(LA(3)^2+4*LA(3)*LS*LL0))]
% ME
[0.5*(lA(4)-sqrt(lA(4)^2+4*LA(4)*LL0*LS)), 0.5*(LA(4)-sqrt(LA(4)^2+4*lA(4)*lS*ll0)),lA(4),0.5*(LA(4)+sqrt(LA(4)^2+4*LA(4)*LS*LL0))]
% R
[0.5*(lA(5)-sqrt(lA(5)^2+4*LA(5)*LL0*LS)), 0.5*(LA(5)-sqrt(LA(5)^2+4*lA(5)*lS*ll0)),lA(5),0.5*(LA(5)+sqrt(LA(5)^2+4*LA(5)*LS*LL0))]

% L3
[0.5*(lA(1)-sqrt(lA(1)^2+4*LA(1)*LL3*LS)), 0.5*(LA(1)-sqrt(LA(1)^2+4*lA(1)*ll3)),lA(1),0.5*(LA(1)+sqrt(LA(1)^2+4*LA(1)*LS*LL3))]
% Rblock
[0.5*(lA(2)-sqrt(lA(2)^2+4*LA(2)*LL3*LS)), 0.5*(LA(2)-sqrt(LA(2)^2+4*lA(2)*ll3)),lA(2),0.5*(LA(2)+sqrt(LA(2)^2+4*LA(2)*LS*LL3))]
% RR
[0.5*(lA(3)-sqrt(lA(3)^2+4*LA(3)*LL3*LS)), 0.5*(LA(3)-sqrt(LA(3)^2+4*lA(3)*ll3)),lA(3),0.5*(LA(3)+sqrt(LA(3)^2+4*LA(3)*LS*LL3))]
% ME
[0.5*(lA(4)-sqrt(lA(4)^2+4*LA(4)*LL3*LS)), 0.5*(LA(4)-sqrt(LA(4)^2+4*lA(4)*ll3)),lA(4),0.5*(LA(4)+sqrt(LA(4)^2+4*LA(4)*LS*LL3))]
% R
[0.5*(lA(5)-sqrt(lA(5)^2+4*LA(5)*LL3*LS)), 0.5*(LA(5)-sqrt(LA(5)^2+4*lA(5)*ll3)),lA(5),0.5*(LA(5)+sqrt(LA(5)^2+4*LA(5)*LS*LL3))]


%% Compute eigenvalues
A = speye((2*N+p)*(k));
A(N*(k)+1:(N+p)*(k),N*(k)+1:(N+p)*(k)) = kron(speye(k),R);
A((N+p)*(k)+1:end,1:N*k) = L';
A(1:N*k,(N+p)*(k)+1:end) = L;
A((N+p)*(k)+1:end,N*k+1:(N+p)*k) = kron(speye(k),Amat.H');
A(N*k+1:(N+p)*k,(N+p)*(k)+1:end) = kron(speye(k),Amat.H);
A((N+p)*(k)+1:end,(N+p)*(k)+1:end) = 0*speye(N*k);

S0 = speye(N*k);

P01 = speye((2*N+p)*(k)); P01(N*(k)+1:(N+p)*(k),N*(k)+1:(N+p)*(k)) = kron(speye(k),Rdiag);
P02 = speye((2*N+p)*(k)); P02(N*(k)+1:(N+p)*(k),N*(k)+1:(N+p)*(k)) = kron(speye(k),Rblock);
P03 = speye((2*N+p)*(k)); P03(N*(k)+1:(N+p)*(k),N*(k)+1:(N+p)*(k)) = kron(speye(k),RRR);
P04 = speye((2*N+p)*(k)); P04(N*(k)+1:(N+p)*(k),N*(k)+1:(N+p)*(k)) = kron(speye(k),RME);
P05 = speye((2*N+p)*(k)); P05(N*(k)+1:(N+p)*(k),N*(k)+1:(N+p)*(k)) = kron(speye(k),R);
%%
P01eig = sort(real(eig(full(P01\A))));
P02eig = sort(real(eig(full(P02\A))));
P03eig = sort(real(eig(full(P03\A))));
P04eig = sort(real(eig(full(P04\A))));
P05eig = sort(real(eig(full(P05\A))));
%%
S1 = inv(L3*L3');

P01((N+p)*(k)+1:end,(N+p)*(k)+1:end) = S1;
P02((N+p)*(k)+1:end,(N+p)*(k)+1:end) = S1;
P03((N+p)*(k)+1:end,(N+p)*(k)+1:end) = S1;
P04((N+p)*(k)+1:end,(N+p)*(k)+1:end) = S1;
P05((N+p)*(k)+1:end,(N+p)*(k)+1:end) = S1;

P31eig = sort(real(eig(full(P01\A))));
P32eig = sort(real(eig(full(P02\A))));
P33eig = sort(real(eig(full(P03\A))));
P34eig = sort(real(eig(full(P04\A))));
P35eig = sort(real(eig(full(P05\A))));