%this routine plots the value function for different values of kappa using
%a simple but inefficient loop. I haven't been able to plot all the value
%functions on the same plot
%% Housekeeping
clear all
close all
%% Defining variables
global r R b kappa z sigma lambda1 lambda2 param_V param_M param_n
r=0.005;                                                               % Discount rate of the agent
R=0.05;                                                                % Opportunity (flow) cost
b=0.05;                                                                
z=0.2;
sigma=0.05;
kappa_vec=[0.15 0.151 0.1511 0.1512 0.1513 0.1514 0.1515 0.1516 0.1517 0.151725 0.15175 0.151775 0.1518 0.1519 0.152 0.153 0.154 0.155 0.156];
for i=1:length(kappa_vec)
kappa=kappa_vec(i);
lambda1=sqrt(2*(r+kappa)/sigma^2);
lambda2=-sqrt(2*(r+kappa)/sigma^2);
param_V=[((kappa*b)/r)   R/r    0   R/r    b/r];                       % Parameters given to fun_findV to find V
param_M=[0                1     0   1       0];                        % Parameters given to fun_findV to find M
param_n=[kappa            0     0   0       1];                        % Parameters given to fun_findV to find N

% Initial conditions for the inaction set: standard Miller-Orr.
mstar_init=((3*b*sigma^2)/(4*R))^(1/3);
mbar_init=3*mstar_init;
policy_init=[mstar_init,mbar_init];

% Constrained minimization to find optimal [mstar,mbar]
LB=[0,0];                                                                %Impose nonnegativity of trigger and target
options = optimset('Display','iter','TolX',1e-10,'FunValCheck','on',...   %Optimization options 
    'TolCon',1e-5);   %Display each iteration
disp('Constrained minimization in progress.');
[policy_opt,fval,EXITFLAG,OUTPUT]=fmincon(@(policy) ...                  %Optimization: constrained minimization
    fun_vstar(policy,param_V), policy_init,[1,-1],0,[],[],LB,[],...      %subject to nonnegativity constraint and
    @constraint,options );                                               %mbar>mstar.

% Finding the optimal policy [mstar_opt,mbar_opt]
mstar_opt=policy_opt(1);
mbar_opt=policy_opt(2);
fprintf('Optimization ended with mstar=%g, mbar=%g and Vstar=%g. \n',mstar_opt,mbar_opt,fval);
% Plotting the value function found
n_eval=1000;
m_vec=0:mbar_opt/(n_eval-1):mbar_opt;
[Vstar,B001,B002,~,~,A0,D0,A,D,B01,B02,B1,B2]=fun_vstar(policy_opt,param_V);
[valuef]=fun_F(B001,B002,A0,D0,A,D,B01,B02,B1,B2,m_vec);
value_function=valuef;
figure(i)
plot(m_vec,valuef,'LineWidth',2);
grid on
hline(Vstar,'--r')
vline(0,'--r');
vline(mstar_opt,'r'); 
vline(mbar_opt,'--r'); 
xlabel('$m$','Interpreter','latex')
ylabel('$V(m,m^{\ast},\bar{m})$','Interpreter','latex')
str = sprintf('Value function for $\\kappa=%g, z=%g, r=%g, R=%g, b=%g, \\sigma=%g$', kappa,z,r,R,b,sigma);
title(str,'Interpreter','latex')

figure(length(kappa_vec)+1)
plot(m_vec,valuef,'LineWidth',2);
hold on
end
hold off
for i=1:1:length(kappa_vec)+1
  saveas(figure(i),fullfile('C:\Users\Elia\Dropbox\RoME\MScThesis\Plots\PlotsModel',['figureValueFKappa' num2str(i) '.png']))
end
