function [path]=simulate(policy,T)
global sigma kappa z 
mstar=policy(1);
mbar=policy(2);
D=1/365;
path=zeros(round(T/D),1);
path(1)=(mbar-mstar)*rand+mstar;
for t=1:round(T/D)
    a1=randn;
    a2=rand;
    if a2<kappa
        if path(t)>z
            path(t+1)=path(t)-z;
        else
            path(t+1)=mstar;
        end
    else
        if path(t)+sigma*a1*sqrt(D)>=mbar || path(t)+sigma*a1*sqrt(D)<=0
            path(t+1)=mstar;
        else
            path(t+1)=path(t)+sigma*a1*sqrt(D);
        end
    end
end
