disp('Testing for Proposition 2 ...');

% Test results from Proposition 2
% Assumption: Proposition 2 requires norm(AA') \le 1
% Parameters
N = 20;
nsubw = 50;
k = 3;

Eval = rand(N,1);
V = rand(N);
A = 1/100*V*diag(Eval)*V';

C = diag(-1*ones(nsubw,1),-1);
L = kron(eye(nsubw+1),eye(N))+kron(C,A);

% Assemble L_M
LMexact = Linvfun(A,nsubw,k);
Lprod2 = (L*LMexact)'*(L*LMexact);

% Define A_1, A_2, A_3
A1 = 0*Lprod2;
Ltemp = Lprod2-eye(N*(nsubw+1));
fin =floor(nsubw/k); 2;
for inc = 1:fin
    (inc-1)*k*N+1
    A1((inc-1)*k*N+1:inc*k*N,(inc-1)*k*N+1:inc*k*N) = ...
        Ltemp((inc-1)*k*N+1:inc*k*N,(inc-1)*k*N+1:inc*k*N);
end
A23 = Ltemp-A1;
A3 = 0*A23;
for inc = 1:2:fin-2
    A3(inc*k*N+1:((inc+1)*k+1)*N,inc*k*N+1:((inc+1)*k+1)*N) = ...
        A23(inc*k*N+1:((inc+1)*k+1)*N,inc*k*N+1:((inc+1)*k+1)*N);
    inc
end
A2 = A23-A3;

% Claim in Proposition 2: eigenvalues of blocks are equivalent to those of
% M^6+M^4+M^2 for k=3, and up to M^2k for other choices of k
figure(3);clf
[V1,E1] = eig(A1(1:k*N,1:k*N));
[V3,E3] = eig(A1(k*N+1:2*k*N,k*N+1:2*k*N));
[V5,E5] = eig(A1(2*k*N+1:3*k*N,2*k*N+1:3*k*N));
sumval = 0*A;
for inc = 1:k
    sumval = sumval+A^(2*inc);
end

[V2,E2] = eig(sumval);
semilogy(diag(E1));
hold on
semilogy(diag(E3));
semilogy(diag(E5));
semilogy([N*(k-1)+1:N*k],sort(diag(E2)));

% Off-diagonal blocks
[V0,E0] = eig(A2(1:(k+1)*N,1:(k+1)*N));
E0 = sort(diag(E0));
% Plot eigenvalues
semilogy([N*(k-1)+1:N*k],E0(N*k+1:N*(k+1)).^2);
% Figure should have all the lines on top of one another (ignore negative
% eigenvalues, though they are also correct)


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