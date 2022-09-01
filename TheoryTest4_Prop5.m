disp('Testing for Proposition 5 ...');

% Test results from Proposition 5
% Assumptions: M must be symmetric with norm \le 1; nsubw = 3,4,5; k = 3
N = 20;
nsubw = 3;
k = 3;

Eval = rand(N,1);
V = rand(N);
A = 1/100*V*diag(Eval)*V';
A = 0.5*(A+A');

% Does A satisfy the characteristic equation?
C = diag(-1*ones(nsubw,1),-1);
L = kron(eye(nsubw+1),eye(N))+kron(C,A);

% Assemble L_M
LMexact = Linvfun(A,nsubw,k);
Lprod2 = (L*LMexact)'*(L*LMexact);

muvec = eig(A);
evec = eig(Lprod2);
nuvec = nufun(muvec);

% Test for extreme eigenvalues (use to produce bounds)
t = 20;
umat = Amatmu(muvec(t));
Prop5_test1 = eig(umat)'+1
Prop5_test2 = nufun(muvec(t))
figure(4); clf
plot(sort(real(evec)),'b-')
hold on
plot(nuvec(1:N),'gx')
plot([N*nsubw+1:N*(nsubw+1)],nuvec(N+1:2*N),'rx')
% We should obtain coincidence of non-zero eigenvalues


function mat = Amatmu(mu)
mat = [mu^6,mu^5,mu^4,-mu^3;
    mu^5,mu^4,mu^3,-mu^2;
    mu^4,mu^3,mu^2,-mu;
    -mu^3,-mu^2,-mu,0];
end

function nu = nufun(mu)
musum = mu.^2+mu.^4+mu.^6;

nup = 1+0.5*(musum+sqrt(4*musum+musum.^2));
num = 1+0.5*(musum-sqrt(4*musum+musum.^2));
nu = sort([num',nup']);
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