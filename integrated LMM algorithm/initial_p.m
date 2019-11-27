function param=initial_p(X,n) %%input X??,n?????
param = struct();

% x1=max(X)
% x2=min(X)
% x3=(x2-x1)*rand(1,n)
sigmazero=zeros(1,n)
muzero=zeros(1,n)
wzero=zeros(1,n)
x5=rand(1,n)
x6=sum(x5)
x7=x5./x6
% x8=1/n
for i1=1:n
    x4=randsrc(1,20,X')
    x1=std(log(x4))
    x2=mean(log(x4))
    sigmazero(1,i1)=x1
    muzero(1,i1)=x2
    wzero(1,i1)=1/n
end
param.mu=muzero
param.sigma=sigmazero
param.weight=wzero



end