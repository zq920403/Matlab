
tic
% first run code in all integrate LMM 
thr=0.015 %set the terminal threshold 
n=3   % set the initial k-value 
X=X8  % X is the travel time dataset
subs=X
pa=initial_p(subs,n)% initialize LMM
[subs,distance,para]=estep(subs,pa,n,1) % run the e-step in EM algorithm
iter=0
while abs(distance(1)-distance(2))>thr & iter<500 % two terminanl condition, one is threshold, the other is max iteration
    pa1=mstep(subs,n)
    [subs,distance,para]=estep(subs,pa1,n,0,distance)
    iter=iter+1
end
%% calculate K-S test for the fitted LMM above 
    nbin=60 
    h=histogram(X,nbin) %% plot histogram of dataset 
    sp=h.BinLimits(1):h.BinWidth:h.BinLimits(2)	
    if length(sp)~=nbin+1
        sp(nbin+1)=h.BinLimits(2)
    end
    sp=sp'
 for i2=1:n
            q(i2)=makedist('lognormal','mu',para.mu(i2),'sigma',para.sigma(i2))
            weight(i2)=para.weight(i2)
            X1(:,i2+2)=weight(i2)*cdf(q(i2),X)
            sp(:,i2+2)=weight(i2)*cdf(q(i2),sp(:,1))
%             qcelln{i2}=@(t)(weight*cdf(q,t))
 end
  h1=sp(:,2) 
  %% calculate AIC 
    for i=1:nbin
        h2(i)=h1(i+1)-h1(i)  
        h3(i)=h.BinCounts(i)
        h4(i)=(h2(i)*length(X)-h3(i))^2
    end
    sse1=sum(h4,2)
   
    s=length(X)
    L=(-(s/2)*log(2*pi))-(s/2)*log(sse1/s)-(s/2)
    aic1=(2*(3*n-1)-2*L)
    %% run K-S test
    y=[X,X1(:,2)]
    [hh,p,sta,cv]=kstest(X,'CDF',y)
    
    %% run the bisection method to calculate the BI
    aver=mean(X)
    datasubs=X
    r1cv=0.95;
    iter1=0;
    datamin=min(datasubs);datamax=max(datasubs);
    r1=@(t)(weight(1)*cdf(q(1),t)+weight(2)*cdf(q(2),t)+weight(3)*cdf(q(3),t))
    while (r1(datamin)<r1cv && r1(datamax)>r1cv &&(r1(datamax)-r1cv)>0.00001)
        iter1=iter1+1;
     midvalue=r1((datamin+datamax)/2);
         if (midvalue<r1cv)
             datamin=(datamin+datamax)/2;
         else 
            datamax=(datamin+datamax)/2;
         end
    
    end
    datamax=double(datamax) 
    BI=(datamax-aver)/aver
 toc