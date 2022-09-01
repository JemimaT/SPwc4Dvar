% number of evaluations for R, B/Q, Rinv, B/Qinv, M in forward and M in
% precond for Lorenz 96 problem

clear
load('14thJanExp.mat','itrand*') %loads iteration info for Lorenz96 experiment


N = 5;
% k =1;
for inc = 1:16
    [R,~,D,Dinv,M,Mpre] = Matvaldiagout(itranddiag(inc,1),15,inc-1,0);%PD.itdiagk1(:,1),5,1,0); % bad R precond
    Rvec(inc) = R;
    Dvec(inc) = D;
    Dinvvec(inc) = Dinv;
    if inc >1
        Mvec(inc) = M+Mpre;
    else 
        Mvec(inc) = M;
    end
    [R,Rinv,D,Dinv,M,Mpre] = Matvaldiagout(itranddiag(inc,2),15,inc-1,1);
    Rvec2(inc) = R;
    Rinvvec2(inc) = Rinv;
    Dvec2(inc) = D;
    Dinvvec2(inc) = Dinv;
    if inc>1
        Mvec2(inc) = M+Mpre;
    else 
        Mvec2(inc) = M;
    end
    % ic
    [R,~,D,M,Mpre] = MatvalICout(itrandic(inc,1),15,inc-1,0);
    RvecIC(inc) = R;
    DvecIC(inc) = D;
    if inc>1
        MvecIC(inc) = M+Mpre;
    else
        MvecIC(inc) = M;
    end
    [R,Rinv,D,M,Mpre] = MatvalICout(itrandic(inc,2),15,inc-1,1);
    Rvec2IC(inc) = R;
    Rinvvec2IC(inc) = Rinv;
    Dvec2IC(inc) = D;
    if inc>1
        Mvec2IC(inc) = M+Mpre;
    else
        Mvec2IC(inc) = M;
    end
end
VarNames = {'R',  'D', 'Dinv','M','R2','Rinv2','D2','Dinv2','M2'};
T1 = table(Rvec',Dvec',Dinvvec',Mvec',Rvec2',Rinvvec2',Dvec2',Dinvvec2',Mvec2', 'VariableNames',VarNames)
%%
VarNames = {'R',  'D','M','R2','Rinv2','D2','M2'};
T2 = table(RvecIC',DvecIC',MvecIC',Rvec2IC',Rinvvec2IC',Dvec2IC',Mvec2IC', 'VariableNames',VarNames)



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