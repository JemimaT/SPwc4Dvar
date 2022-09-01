disp('Testing for Proposition 1 ...');

% Test result and intermediary steps of Proposition 1
% N.B.: we do not have any assumptions on the structure of M
load icholcompN250.mat

% Parameter selection
nsubw = 13; % number of blocks
k = 8; % skipping parameter in the preconditioner

M = Mtlen+eye(250); N = length(M);

% Assemble L
C = diag(-1*ones(nsubw,1),-1);
L = kron(eye(nsubw+1),eye(N))+kron(C,M);

% Assemble L_M
LMexact = Linvfun(M,nsubw,k);
Lprod2 = (L*LMexact)'*(L*LMexact);

[evec,eval] = eig(Lprod2-eye(N*(nsubw+1)));
eval = diag(eval);
% Test the result of Proposition 1
Prop1_test1 = sum(abs(eval)<1e-10);
r = nsubw+1-2*floor(nsubw/k);
Prop1_test1 = [Prop1_test1 r*N]
% the two printed values should be the same, apart from numerical round-off

% Show ansatz of eigenvectors of A(M) in Prop 1 relates to zero eigenvector
diffvec = zeros([r*N,1]);
% Zero eigenvalues from first block
for i = 1:N
   diffvec(i) = ...
       sum(abs(testvecfun(M,nsubw,i)-Lprod2*testvecfun(M,nsubw,i)));
end
ind = N+1;
for inc = 1:floor(nsubw/k) 
    for i = (k*(inc-1)+1)*N+1:(inc*(k)-1)*N
        diffvec(ind) = sum(abs((Lprod2-eye(N*(nsubw+1)))* ...
            testvecfunmult(M,nsubw,i,1+floor((i-1)/N))));
        ind = ind+1;
    end
end

% Check the additional unit eigenvectors
In = 0*Lprod2;
In(N*nsubw+1:N*(nsubw+1),N*nsubw+1:N*(nsubw+1)) = eye(N);
for inc = 1:nsubw-k*floor(nsubw/k)
    for i = 1:N
        diffvec(ind) = sum(abs((Lprod2-eye(N*(nsubw+1)))*In(:,i)));
        ind = ind+1;
    end
end

% Check number of eigenvalues agrees with result in proof; want close to 0
Prop1_test2 = r*N-ind+1 

% Check maximum sum of abs(sum(Av)) (should return zero due to zero vector)
Prop1_test3 = max(diffvec) % want close to zero


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

function evec = testvecfun(A,nsubw,ei)
% Eigenvectors for first block
N = length(A);
evec = zeros([N*(nsubw+1),1]);
evec(ei) = 1;
remi = ei;
vec2 = zeros([N,1]); vec2(remi) = 1;
Ai = A*vec2;
evec(1+N:2*N) = -Ai;
end

function evec = testvecfunmult(A,nsubw,ei,block2)
% Eigenvectors for remaining block
% ei is entry for canonical eigenvector in v_t
% block2 tells us where to read in value of -Av_t
N = length(A);
evec = zeros([N*(nsubw+1),1]);
evec(ei) = 1;
remi = ei-N*floor(ei/N);
if remi ==0
    remi = N;
end
vec2 = zeros([N,1]); vec2(remi) = 1;
Ai = A*vec2;
evec(1+(block2)*N:(block2+1)*N) = -Ai;
end