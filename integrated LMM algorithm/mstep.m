function M=mstep(X,n)
    S=length(X)
    for i3=1:n
        sumRIK(1,i3)=sum(X(:,i3+2),1)
        w1(1,i3)= sumRIK(1,i3)./S  %%w1 is new weight vector  
        munew(1,i3)=sum(log(X(:,1)).*X(:,i3+2),1)./sumRIK(1,i3)
        sigma2new(1,i3)= (sum(X(:,i3+2).*((log(X(:,1))- munew(1,i3)).^2),1))./sumRIK(1,i3)
    end
    M=struct();
    M.mu=munew
    M.sigma=sqrt(sigma2new)
    M.weight=w1
   
end