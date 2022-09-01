disp('Testing for Lemma 1 ...');

N = 15; % N+1 blocks
np = 10; % number of points per block
k = 4; % skip every k blocks

Md = rand(N*np,np); % consists of lower triangular blocks of L

L = eye((N+1)*np);
for i = 1:N
   L(i*np+1:(i+1)*np,(i-1)*np+1:i*np) = -Md((i-1)*np+1:i*np,:); 
end

LM = L;
for i = 1:floor(N/k)
    LM(i*k*np+1:i*k*np+np,i*k*np-np+1:i*k*np) = 0;
end

LMinv = eye((N+1)*np);
for i = 1:N+1
    for j = 1:N+1
        if (i-j >= 1) && (i-j <= mod(i-1,k))
            for m = 1:i-j
                if m==1
                    LMinv((i-1)*np+1:i*np,(j-1)*np+1:j*np) = eye(np);
                    LMinv((i-1)*np+1:i*np,(j-1)*np+1:j*np) = ...
                        LMinv((i-1)*np+1:i*np,(j-1)*np+1:j*np)* ...
                        Md((i-m-1)*np+1:(i-m)*np,:);                    
                else
                    LMinv((i-1)*np+1:i*np,(j-1)*np+1:j*np) = ...
                        LMinv((i-1)*np+1:i*np,(j-1)*np+1:j*np)* ...
                        Md((i-m-1)*np+1:(i-m)*np,:);
                end
            end
        end
    end
end

Lem1_test1 = norm(LMinv-inv(LM),inf)/norm(LMinv,inf) % want close to zero
Lem1_test2 = norm(LM*LMinv-eye((N+1)*np),inf) % want close to zero
figure(1), surf(LM*LMinv-eye((N+1)*np)) % want close to identity matrix

%%
disp('Testing for Lemma 2 ...');

A = zeros((N+1)*np,(N+1)*np);
for n = 1:floor(N/k) % second term of A
    i = n*k+1;
    for j = (n-1)*k+1:n*k
        A((i-1)*np+1:i*np,(j-1)*np+1:j*np) = -eye(np);
        for t = j:n*k
            A((i-1)*np+1:i*np,(j-1)*np+1:j*np) = ...
                A((i-1)*np+1:i*np,(j-1)*np+1:j*np)* ...
                Md((n*k-t+j-1)*np+1:(n*k-t+j)*np,:);
        end
    end
end
A = A+A'; % to take account of third term as well

A1stTerm = zeros((N+1)*np,(N+1)*np);
for n = 1:floor(N/k) % first term of A
    for i = (n-1)*k+1:n*k
        for j = (n-1)*k+1:n*k
            A1stTerm((i-1)*np+1:i*np,(j-1)*np+1:j*np) = eye(np);
            for t = i:n*k
                A1stTerm((i-1)*np+1:i*np,(j-1)*np+1:j*np) = ...
                    A1stTerm((i-1)*np+1:i*np,(j-1)*np+1:j*np)* ...
                    Md((t-1)*np+1:t*np,:)';
            end
            for q = j:n*k
                A1stTerm((i-1)*np+1:i*np,(j-1)*np+1:j*np) = ...
                    A1stTerm((i-1)*np+1:i*np,(j-1)*np+1:j*np)* ...
                    Md((n*k-q+j-1)*np+1:(n*k-q+j)*np,:);
            end
        end
    end
end

A = eye((N+1)*np)+A1stTerm+A;

Lem2_test = norm(LMinv'*(L'*L)*LMinv-A,inf)/norm(A,inf) % want close to 0
figure(2), surf(LMinv'*(L'*L)*LMinv-A) % want close to zero matrix