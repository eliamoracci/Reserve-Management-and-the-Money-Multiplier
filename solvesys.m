function x=solvesys(par)
global lambda1 lambda2 z
A=cell2mat(par(1));
D=cell2mat(par(2));
B01=cell2mat(par(3));
B02=cell2mat(par(4));
B1=cell2mat(par(5));
 B2=cell2mat(par(6));
j=cell2mat(par(7));
sum1=0;
for i=1:j %#ok<*IJCL>
    sum1=sum1+B1(j,i)*exp(lambda1*z)*(z^i)+B2(j,i)*exp(lambda2*z)*(z^i);
end
sum2=0;
for i=1:j %#ok<*IJCL>
    sum2=sum2+B1(j,i)*exp(lambda1*z)*(lambda1*(z^i)+i*(z^(i-1)))+B2(j,i)*exp(lambda2*z)*(lambda2*z^i+i*(z^(i-1)));
end
matA=[1,1;lambda1,lambda2];
vecb=[A(j)+D(j)*z+B01(j)*exp(lambda1*z)+B02(j)*exp(lambda2*z)+sum1-A(j+1);...
    D(j)+B01(j)*exp(lambda1*z)*lambda1+B02(j)*exp(lambda2*z)*lambda2+sum2-D(j+1)-B1(j+1,1)-B2(j+1,1)];
x=linsolve(matA,vecb);
end