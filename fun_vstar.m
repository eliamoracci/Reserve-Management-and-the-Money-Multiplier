%This function solves the level of the value function at mstar,
%V(mstar)=Vstar, for a given value mstar and mbar (vector policy). 
%To solve for Vstar, one should guess Vstar and a parameter of the solution, 
%B001 (see Appendix of %Alvarez-Lippi paper and the description of the 
%function fun_findV). 

%The function fsolve finds the value of Vstar and B001 so that Vstar is equal
%to imp1_Vstar and imp2_Vstar, the two values implied by the shape of the value
%function derived from initial parameters Vstar and B001.


function [Vstar,B001,B002,J,Jstar,A0,D0,A,D,B01,B02,B1,B2] = fun_vstar(policy,param)
global kappa b r R sigma 
  
% Initial conditions for Vstar and B001
Vstar_init  = (R*(4/3)*((3*b*sigma^2)/(4*R))^(1/3)+b*(sigma^2/2)*((3*b*sigma^2)/(4*R))^(-2/3))/r+24;
B001_init=(Vstar_init+b-((kappa*(Vstar_init+b))/(r+kappa)))/25;
input_init=[Vstar_init,B001_init];
%Solving for Vstar and B001
options=optimoptions('fsolve','Display','off','StepTolerance',1e-30,'FunValCheck',...
    'on','FunctionTolerance',1e-6,'Algorithm','trust-region-dogleg');
[sol,FVAL_V,EXITFLAG,OUTPUT] = fsolve(@(input) fun_findV(input,policy,param), input_init, options);
%Getting Vstar, the output.
Vstar=sol(1);
[~,B001,B002,J,Jstar,A0,D0,A,D,B01,B02,B1,B2]=fun_findV(sol,policy,param);
end
