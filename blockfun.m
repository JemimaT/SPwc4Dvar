% Partition by testing off diagonal blocks

function [Rapprox,splitting,nproc] = blockfun(R,pvec,tol,maxsize,processor)
%addpath '/home/jtabeart/Documents/EdPDRA/DA_Precond/DA_Precond'

% want to make a function
% inputs matrix, block sizes/pvec, number of processors, tolerance for
% allocation, max size per processor

p = length(R);

%processor = 3;
%maxsize = 3;

% define test R (remove when do loop)
pstart = cumsum(pvec)-pvec+1; %where blocks start

% go through off diagonal blocks
plen = length(pvec); % how many blocks
pstart(plen+1) = pvec(plen)+pstart(plen)-1;
normvec = zeros([1,plen-1]);
for inc = 1:plen-1 % compute norm of each off off-diag blocks
    normvec(inc) = norm(R(pstart(inc):pstart(inc+1)-1,pstart(inc+1):pstart(inc+2)-1),'fro')/sqrt(pvec(inc)*pvec(inc+1));
end
normvec
[v1,v2] = sort(normvec); % sort norms in increasing order
% want to keep norm that is bigger than tolerance

% Check if we have correct number of processors
% decide where to start in v1
%divvec = zeros([1,nnz(v1<tol)]);

divveczero = sort(v2(v1<tol)); % beginning and end of zero blocks
%%
Rapprox = sparse(p,p);
%divloc = find(divvec);
%divvec(1) = 1; divvec(end) = plen;
% if max(diff(divveckeep))>10
%     
%     for inc2 = 1:plen-1
%         Rapprox(pstart(v2(inc2)):pstart(v2(inc2)+1)-1,pstart(v2(inc2)):pstart(v2(inc2)+1)-1)  = R(pstart(v2(inc2)):pstart(v2(inc2)+1)-1,pstart(v2(inc2)):pstart(v2(inc2)+1)-1) ;
%         if v1(inc2)>tol
%             Rapprox(pstart(v2(inc2)):pstart(v2(inc2)+1)-1,pstart(v2(inc2)+1):pstart(v2(inc2)+2)-1) =R(pstart(v2(inc2)):pstart(v2(inc2)+1)-1,pstart(v2(inc2)+1):pstart(v2(inc2)+2)-1);
%             Rapprox(pstart(v2(inc2)+1):pstart(v2(inc2)+2)-1,pstart(v2(inc2)):pstart(v2(inc2)+1)-1) = R(pstart(v2(inc2)+1):pstart(v2(inc2)+2)-1,pstart(v2(inc2)):pstart(v2(inc2)+1)-1);
%             
%         end
%     end
%     Rapprox(pstart(plen):pstart(plen+1),pstart(plen):pstart(plen+1)) = R(pstart(plen):pstart(plen+1),pstart(plen):pstart(plen+1));
% else
Rapprox(pstart(1):pstart(divveczero(1)+1)-1,pstart(1):pstart(divveczero(1)+1)-1)= R(pstart(1):pstart(divveczero(1)+1)-1,pstart(1):pstart(divveczero(1)+1)-1);
splitting = zeros([1,length(divveczero)+1]);
splitting(1) = pstart(divveczero(1)+1)-1;
for inc2 = 1:length(divveczero)-1
    
    
    Rapprox(pstart(divveczero(inc2)+1):pstart(divveczero(inc2+1)+1)-1,pstart(divveczero(inc2)+1):pstart(divveczero(inc2+1)+1)-1) = R(pstart(divveczero(inc2)+1):pstart(divveczero(inc2+1)+1)-1,pstart(divveczero(inc2)+1):pstart(divveczero(inc2+1)+1)-1);
    splitting(inc2+1) = pstart(divveczero(inc2+1)+1)-(pstart(divveczero(inc2)+1));
end
Rapprox(pstart(divveczero(end)+1):end,pstart(divveczero(end)+1):end) = R(pstart(divveczero(end)+1):end,pstart(divveczero(end)+1):end);
splitting(end)= p-pstart(divveczero(end)+1)+1;
pstart2 = ones([1,length(splitting)+1]);
pstart2(2:end) = cumsum(splitting);

% Check if we have small enough blocks
while max(splitting)>maxsize
    %split
    [w1,w2] = max(splitting);
    % divide this block in two
    s1 = floor(w1/2);s2 = w1-floor(w1/2);
    % set off-diag  blocks to zero
    Rapprox(pstart2(w2):pstart2(w2)+s1,pstart2(w2)+s1+1:pstart2(w2)+w1) = zeros(s1+1,s2); 
    Rapprox(pstart2(w2)+s1+1:pstart2(w2)+w1,pstart2(w2):pstart2(w2)+s1) = zeros(s2,s1+1); 
    
    % update splitting
    splittingtemp = splitting;
    splitting = zeros([1,1+length(splittingtemp)]);
    splitting(1:w2-1)= splittingtemp(1:w2-1);
    splitting(w2) = floor(w1/2);
    splitting(w2+1) = w1-floor(w1/2);
    splitting(w2+2:end) = splittingtemp(w2+1:end);
end %do nothing
pstart2 = ones([1,length(splitting)+1]);
pstart2(2:end) = cumsum(splitting);
nproc = length(splitting);
while nproc > processor
    comb = splitting(1:end-1)+splitting(2:end);
    [~,m] = min(comb); %which two to combine % m and m+1th
    
    % combine smallest blocks
    Rapprox(pstart2(m)+1:pstart2(m+2),pstart2(m)+1:pstart2(m+2)) = R(pstart2(m)+1:pstart2(m+2),pstart2(m)+1:pstart2(m+2));
    % update splitting
    splittingtemp = splitting;
    splitting = zeros([1,length(splittingtemp)-1]);
    splitting(1:m-1) = splittingtemp(1:m-1);
    splitting(m) = splittingtemp(m)+splittingtemp(m+1);
    if m<length(splittingtemp)-1
        splitting(m+1:end) = splittingtemp(m+2:end);
    end
    nproc = nproc-1;
    
end

end
