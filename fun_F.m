%This function computes the value of the function V at points specified by
%the vector m_vec. The output F_m is a vector which contains V(m) for
%m=[mvec]. The inputs are the coefficients of the value function and the
%vector m_vec.
function [F_m]=fun_F(B001,B002,A0,D0,A,D,B01,B02,B1,B2,m_vec)
global lambda1 lambda2 z 
F_m=zeros(length(m_vec),1);
for i=1:length(m_vec)
    j=round(m_vec(i)/z);
    if j>m_vec(i)/z
        j=j-1;
    end
    if j==0
        F_m(i)=A0+D0*m_vec(i)+B001*exp(lambda1*m_vec(i))+B002*exp(lambda2*m_vec(i));
    else
        F_m(i)=A(j)+D(j)*(m_vec(i)-z*j)+B01(j)*exp(lambda1*(m_vec(i)-z*j))+B02(j)*exp(lambda2*(m_vec(i)-z*j));
        for h=1:j
            F_m(i)=F_m(i)+B1(j,h)*exp(lambda1*(m_vec(i)-z*j))*((m_vec(i)-z*j)^h)+...
                +B2(j,h)*exp(lambda2*(m_vec(i)-z*j))*((m_vec(i)-z*j)^h);
        end
    end
end