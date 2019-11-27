function [E,distance1,param1]=estep(X,param,n,Isfirst,distance)
    for i2=1:n
            q=makedist('lognormal','mu',param.mu(i2),'sigma',param.sigma(i2))
            weight=param.weight(i2)
            qcell{i2}=@(t)(weight*pdf(q,t))
            X(:,i2+2)=qcell{i2}(X(:,1))
    end

%     [nrow,ncol]=size(X)
%     for i2=1:n
%         X(:,i2+2)=qcell{i2}(X(:,1))
%     end
    X1=X(:,1:2)
    sum1=sum(X(:,3:2+n),2)
    X2=X./sum1
    X2(:,1:2)=X1
    X=X2
%     for i2=1:n
%         X(:,i2+2)=X(:,i2+2)./sum1
%     end
    [max1,index]=max(X(:,3:2+n),[],2)
    X(:,2)=index
    E=X
    for i2=1:n
        lim(i2,1)=sum(abs(log(X(X(:,2)==i2,1))-param.mu(i2)),1)
    end
    
    if Isfirst==1
       distance1(1,1)=sum(lim,1)
       distance1(1,2)=0
       param1=param
    else
%         distance1=distance
        distance1(1,2)=distance(1,1)
        distance1(1,1)=sum(lim,1)
        param1=param
    end
    
end