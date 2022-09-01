disp('Testing for Proposition 6 ...');

% Test results from Proposition 6
% Assumptions: M must be symmetric with norm \le 1; N = 8,9,10,11; k = 4
N = 20;
nsubw = 10;
k = 4;

Eval = rand(N,1);
V = rand(N);
A = 1/100*V*diag(Eval)*V';
A = 0.5*(A+A');

% Assemble matrices
C = diag(-1*ones(nsubw,1),-1);
L = kron(eye(nsubw+1),eye(N))+kron(C,A);

% Assemble L_M
LMexact = Linvfun(A,nsubw,k);
Lprod2 = (L*LMexact)'*(L*LMexact);

muvec = eig(A);
evec = eig(Lprod2);

t = 20;
umat = Amatmu(muvec(t));
Prop6_test = eig(umat)'+1

figure(5); clf
plot(sort(real(evec)),'b-')
hold on

[mu1,mu2] = Amatsplit(muvec(t));

maxmu1 = eigs(mu1,1,'lr');
maxmu2 = eigs(mu2,1,'lr');
plot([1,length(Lprod2)],[maxmu1+maxmu2+1,maxmu1+maxmu2+1],'r--')
plot([1,length(Lprod2)],[5+sqrt(8),5+sqrt(8)],'k--')


function mat = Amatmu(mu)
mat = [mu^8,mu^7,mu^6,mu^5,-mu^4,0,0,0,0;
    mu^7,mu^6,mu^5,mu^4,-mu^3,0,0,0,0;
    mu^6,mu^5,mu^4,mu^3,-mu^2,0,0,0,0;
    mu^5,mu^4,mu^3,mu^2,-mu^1,0,0,0,0;
    -mu^4,-mu^3,-mu^2,-mu^1,mu^8,mu^7,mu^6,mu^5,-mu^4;
    0,0,0,0,mu^7,mu^6,mu^5,mu^4,-mu^3;
    0,0,0,0,mu^6,mu^5,mu^4,mu^3,-mu^2;
    0,0,0,0,mu^5,mu^4,mu^3,mu^2,-mu;
    0,0,0,0,-mu^4,-mu^3,-mu^2,-mu,0];
end

function [mat1,mat2] = Amatsplit(mu)
mat1 = [mu^8,mu^7,mu^6,mu^5,0,0,0,0,0;
    mu^7,mu^6,mu^5,mu^4,0,0,0,0,0;
    mu^6,mu^5,mu^4,mu^3,0,0,0,0,0;
    mu^5,mu^4,mu^3,mu^2,0,0,0,0,0;
    0,0,0,0,mu^8,0,0,0,0;
    0,0,0,0,0,mu^6,mu^5,mu^4,-mu^3;
    0,0,0,0,0,mu^5,mu^4,mu^3,-mu^2;
    0,0,0,0,0,mu^4,mu^3,mu^2,-mu;
    0,0,0,0,0,-mu^3,-mu^2,-mu,0];

mat2 = [0,0,0,0,-mu^4,0,0,0,0;
    0,0,0,0,-mu^3,0,0,0,0;
    0,0,0,0,-mu^2,0,0,0,0;
    0,0,0,0,-mu^1,0,0,0,0;
    -mu^4,-mu^3,-mu^2,-mu^1,0,mu^7,mu^6,mu^5,-mu^4;
    0,0,0,0,mu^7,0,0,0,0;
    0,0,0,0,mu^6,0,0,0,-0;
    0,0,0,0,mu^5,0,0,0,0;
    0,0,0,0,-mu^4,0,0,0,0];
end

function LMinv = Linvfun(M,nsubw,k)
N = length(M);

floorvar = floor((nsubw+1)/k); % whole blocks excluding the first
Mi = eye(N);
remval = nsubw+1-k*floorvar;
temp = eye(N*k);
for inc1 = 1:k-1
    Mi = Mi*M;
    Ci = diag(ones(k-inc1,1),-inc1);
    temp = temp+kron(Ci,Mi);
end

Mi = eye(N);
remvec = eye(N*remval);
for inc1 = 1:remval-1
    Mi = Mi*M;
    Ci = diag(ones(remval-inc1,1),-inc1);
    remvec = remvec+kron(Ci,Mi);
end
% Assemble matrix
Cfull = ones(floorvar+1,1); Cfull(end) = 0;
temp2 = kron(diag(Cfull),temp);
LMinv = temp2(1:N*(nsubw+1),1:N*(nsubw+1));
LMinv(k*floorvar*N+1:end,k*floorvar*N+1:end) = remvec;
end