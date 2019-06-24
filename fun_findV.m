%This function takes as inputs Vstar, B001, mstar and mbar. For given
%values of these variables this function computes the function V and then
%two implied values for Vstar, derived by the conditions imp1_Vstar=V(mstar) 
%and imp2_Vstar=V(mbar)-b. 

%The two outputs are the distance between imp1_Vstar and imp2_Vstar and the
%distance betwwen imp2_Vstar and the initial Vstar.

%I define the function and all of its inputs.
function [zero,B001,B002,J,Jstar,A0,D0,A,D,B01,B02,B1,B2]=fun_findV(input,policy,param)
global kappa z r sigma lambda1 lambda2
alpha_0 = param(1);
nu_0 = param(2);
alpha   = param(3);
nu  = param(4);
alphastar= param(5);
Vstar=input(1);
B001=input(2);
mstar=policy(1);
mbar=policy(2);
%I derive the number of intervals, based on values for mbar/mstar.
J=round(mbar/z);
if J<mbar/z || J==0
    J=J+1;
end
Jstar=round(mstar/z);
if Jstar>mstar/z
    Jstar=Jstar-1;
end

%Formulas for j=0; B001 is an input so I take it as given
D0=(r/(r+kappa))*nu_0;
A0=((kappa*Vstar)/(r+kappa))+alpha_0*(r/(r+kappa));
B002=Vstar+(r*alphastar)-A0-B001;

%Other coefficients of the solution, stored in matrices and vectors
B1=zeros(J-1,J-1);
B01=zeros(J-1,1);
B2=zeros(J-1,J-1);
B02=zeros(J-1,1);
A=zeros(J-1,1);
D=zeros(J-1,1);

%Here I compute {A_j}, {D_j} and the diagonal {B_{j,j}}
for j=0:J-2
    if j==0  % this if takes care of the case of j=0;
        D(j+1)=((r/(r+kappa))*nu)+(kappa/(r+kappa))*D0;
        A(j+1)=((r*alpha)/(r+kappa))+((kappa*A0-z*kappa*(j+1)*D0)/(r+kappa))+D(j+1)*z*(j+1);
        B1(j+1,1)=-(kappa/(lambda1*(sigma^2)))*B001;
        B2(j+1,1)=-(kappa/(lambda2*(sigma^2)))*B002;
    else
        D(j+1)=((r/(r+kappa))*nu)+(kappa/(r+kappa))*D(j);
        A(j+1)=((r*alpha)/(r+kappa))+((kappa*A(j)-z*kappa*(j+1)*D(j))/(r+kappa))+D(j+1)*z*(j+1);
        B1(j+1,j+1)=-(kappa/(lambda1*(sigma^2)*(j+1)))*B1(j,j);
        B2(j+1,j+1)=-(kappa/(lambda2*(sigma^2)*(j+1)))*B2(j,j);
    end
end


for j=0:J-2
    %Here I impose VM and SP to compute B_{1,0}.
    if j==0
        par={A0,A(1),D0,D(1),B001,B002,B1,B2};
        x=solvesys0(par);
        B01(1)=x(1);
        B02(1)=x(2);
    else
        %Having B_{j+1,j+1} as BCs, I recursively compute B_{j+1,i} until
        %i=2
        for i=j:-1:2
             B1(j+1,i)=-(kappa/(lambda1*(sigma^2)*i))*B1(j,i-1)-(1/(lambda1*2))*B1(j+1,i+1)*(i+1);
             B2(j+1,i)=-(kappa/(lambda2*(sigma^2)*i))*B2(j,i-1)-(1/(lambda2*2))*B2(j+1,i+1)*(i+1);
        end
        %Then, I compute B_{j+1,1}
        B1(j+1,1)=-(kappa/(lambda1*(sigma^2)))*B01(j)-(1/lambda1)*B1(j+1,2);
        B2(j+1,1)=-(kappa/(lambda2*(sigma^2)))*B02(j)-(1/lambda2)*B2(j+1,2);

        %To end up, I compute B_{j+1,0}
        par={A,D,B01,B02,B1,B2,j};
        x=solvesys(par); 
        B01(j+1)=x(1);
        B02(j+1)=x(2);
        
    end
    
end

%After the above loop I have all B1, B2, B01, B02; now I have to compute the implied Vstars

if Jstar==0
    imp1_Vstar = A0 + D0*(mstar) + B001*exp(lambda1*mstar)+B002*exp(lambda2*mstar);
elseif Jstar>0 && Jstar<=J-1
    imp1_Vstar = A(Jstar) + D(Jstar)*(mstar-z*(Jstar)) + B01(Jstar)*exp(lambda1*(mstar-z*(Jstar)))+...
            +B02(Jstar)*exp(lambda2*(mstar-z*(Jstar))) ;
        for i=1:Jstar
        imp1_Vstar= imp1_Vstar + B1(Jstar,i)*exp(lambda1*(mstar-z*(Jstar))) * (mstar-z*(Jstar))^i +...
                +B2(Jstar,i)*exp(lambda2*(mstar-z*(Jstar))) * (mstar-z*(Jstar))^i;
        end
elseif Jstar>J-1
   fprintf('Mistake: Jstar has value %d and J has value %d , because mstar=%d and mbar=%d.\n',Jstar,J,mstar,mbar)
end

if J==1
    imp2_Vstar = A0+D0*(mbar)+B001*exp(lambda1*mbar)+B002*exp(lambda2*mbar)-(r*alphastar);
elseif J>1
    imp2_Vstar =  A(J-1) + D(J-1)*(mbar-z*(J-1)) + B01(J-1,1)*exp(lambda1*(mbar-z*(J-1)))+...
            +B02(J-1,1)*exp(lambda2*(mbar-z*(J-1)))-(r*alphastar);
        for i=1:J-1
        imp2_Vstar= imp2_Vstar + B1(J-1,i)*exp(lambda1*(mbar-z*(J-1))) * (mbar-z*(J-1))^i +...
                +B2(J-1,i)*exp(lambda2*(mbar-z*(J-1))) * (mbar-z*(J-1))^i;
        end
else
    fprintf('Mistake: J has value %d because mbar is %d and z is %d.\n',Jstar,mbar,z)
        
end

% Outputs of the function
zero(1)=(imp1_Vstar-Vstar)/r;
zero(2)=(imp2_Vstar-Vstar)/r;
