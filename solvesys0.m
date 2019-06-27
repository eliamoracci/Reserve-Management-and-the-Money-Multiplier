function x=solvesys0(par)
global z lambda1 lambda2
A0=cell2mat(par(1));
A1=cell2mat(par(2));
D0=cell2mat(par(3));
D1=cell2mat(par(4));
B001=cell2mat(par(5));
B002=cell2mat(par(6));
B1=cell2mat(par(7));
B2=cell2mat(par(8));
matA=[1,1;lambda1,lambda2];
vecb=[A0+D0*z+B001*exp(lambda1*z)+B002*exp(lambda2*z)-A1; D0+B001*exp(lambda1*z)*lambda1+B002*exp(lambda2*z)*lambda2-D1-B1(1,1)-B2(1,1)];
x=linsolve(matA,vecb);
end