function [c,ceq]=constraint(policy)
global z lambda1 lambda2 param_V
[~,B001,B002,J,~,~,D0,~,D,B01,B02,B1,B2]=fun_vstar(policy,param_V);
if J>1
    ceq=D(J-1)+B01(J-1)*exp(lambda1*(policy(2)-z*(J-1)))*(lambda1)+B02(J-1)*exp(lambda2*(policy(2)-z*(J-1)))*(lambda2);
    for i=1:(J-1)
        ceq=ceq+B1(J-1,i)*exp(lambda1*(policy(2)-z*(J-1)))*(lambda1*((policy(2)-z*(J-1))^(i))+i*((policy(2)-z*(J-1))^(i-1)))+...
            +B2(J-1,i)*exp(lambda2*(policy(2)-z*(J-1)))*(lambda2*((policy(2)-z*(J-1))^(i))+i*((policy(2)-z*(J-1))^(i-1)));
    end
else
    ceq=D0+B001*exp(lambda1*(policy(2)-z*(J-1)))*(lambda1)+B002*exp(lambda2*(policy(2)-z*(J-1)))*(lambda2);
end
c=[];