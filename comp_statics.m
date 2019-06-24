%This routine performs comparative statics of M and N with respect to kappa
%and with respect to R. Regarding the CS with respect to R, it is performed
%both with kappa>0 and kappa=0 (standard Miller-Orr model).
%%  Comp statics wrt kappa
clear all
close all
global r R b sigma z kappa lambda1 lambda2 param_V param_M param_n
r=0.005;                                                               % Discount rate of the agent
R=0.05;                                                                % Opportunity (flow) cost
b=0.05;                                                                
sigma=0.05;
z=0.17;
mstar_init=((3*b*sigma^2)/(4*R))^(1/3);
mbar_init=3*mstar_init;
policy_init=[mstar_init,mbar_init];
kappa_vec=0.01:0.01:0.25;
for i=1:length(kappa_vec)
    kappa=kappa_vec(i);
    lambda1=sqrt(2*(r+kappa)/sigma^2);
    lambda2=-sqrt(2*(r+kappa)/sigma^2);
    param_V=[((kappa*b)/r)   R/r    0   R/r    b/r];                       % Parameters given to fun_findV to find V
    param_M=[0                1     0   1       0];                        % Parameters given to fun_findV to find M
    param_n=[kappa            0     0   0       1];     
    fprintf('Calculating statistics and optimal policy for kappa=%g. \n',kappa)
    LB=[0,0]; 
    options = optimset('Display','off','TolX',1e-10,'FunValCheck','on','TolCon',1e-5); %Display each iteration
    [policy_opt,fval,EXITFLAG,OUTPUT]=fmincon(@(policy) fun_vstar(policy,param_V),...
    policy_init,[1,-1],0,[],[],LB,[],@constraint,options );          %Constraints: lower bound and mbar>mstar.
    mstar1(i)=policy_opt(1);
    mbar1(i)=policy_opt(2);
    [V1(i),~,~,~,~,~,~,~,~,~,~,~,~]=fun_vstar(policy_opt,param_V);
    [M1(i),~,~,~,~,~,~,~,~,~,~,~,~]=fun_vstar(policy_opt,param_M);
    [N1(i),~,~,~,~,~,~,~,~,~,~,~,~]=fun_vstar(policy_opt,param_n);
end
figure(1)
plot(kappa_vec,M1,'LineWidth',2,'Color',[0 0.5 1])
hold on
plot(kappa_vec,N1,'LineWidth',2,'Color',[0 0.5 1],'LineStyle',':')
str = sprintf('Aggr. stats for $z=%g, r=%g, R=%g, b=%g, \\sigma=%g$ as $\\kappa$ varies', z,r,R,b,sigma);
title(str,'Interpreter','latex','FontSize',10)
xlabel('$\kappa$','Interpreter','latex')
leg1=legend('$M(\kappa)$','$n(\kappa)$');
set(leg1,'Interpreter','latex','Location','northwest');
figure(2)
plot(kappa_vec,mstar1,'LineWidth',2,'Color',[0 0.6 0],'LineStyle',':')
hold on
plot(kappa_vec,mbar1,'LineWidth',2,'Color',[0 0.6 0])
str = sprintf('Optimal policy for $z=%g, r=%g, R=%g, b=%g, \\sigma=%g$ as $\\kappa$ varies', z,r,R,b,sigma);
title(str,'Interpreter','latex','FontSize',10)
xlabel('$\kappa$','Interpreter','latex')
leg2=legend('$m^{\ast}(\kappa)$','$\bar{m}(\kappa)$');
set(leg2,'Interpreter','latex','Location','northwest');
%% 5) Comp statics wrt R 
r=0.005;                                                               % Discount rate of the agent                                                                % Opportunity (flow) cost
b=0.05;                                                                
sigma=0.05;
z=0.17;
kappa=0.125;
R=0.05;
mstar_init=((3*b*sigma^2)/(4*R))^(1/3);
mbar_init=3*mstar_init;
policy_init=[mstar_init,mbar_init];
R_vec=[0.01:0.005:0.25];
for i=1:length(R_vec)
    R=R_vec(i);
    lambda1=sqrt(2*(r+kappa)/sigma^2);
    lambda2=-sqrt(2*(r+kappa)/sigma^2);
    param_V=[((kappa*b)/r)   R/r    0   R/r    b/r];                       % Parameters given to fun_findV to find V
    param_M=[0                1     0   1       0];                        % Parameters given to fun_findV to find M
    param_n=[kappa            0     0   0       1];     
    fprintf('Calculating statistics and optimal policy for R=%g. \n',R)
    LB=[0,0]; 
    options = optimset('Display','off','TolX',1e-10,'FunValCheck','on','TolCon',1e-5); %Display each iteration
    [policy_opt,fval,EXITFLAG,OUTPUT]=fmincon(@(policy) fun_vstar(policy,param_V),...
    policy_init,[1,-1],0,[],[],LB,[],@constraint,options );          %Constraints: lower bound and mbar>mstar.
    mstar2(i)=policy_opt(1);
    mbar2(i)=policy_opt(2);
    [V2(i),~,~,~,~,~,~,~,~,~,~,~,~]=fun_vstar(policy_opt,param_V);
    [M2(i),~,~,~,~,~,~,~,~,~,~,~,~]=fun_vstar(policy_opt,param_M);
    [N2(i),~,~,~,~,~,~,~,~,~,~,~,~]=fun_vstar(policy_opt,param_n);
end

figure(3)
plot(R_vec,M2,'LineWidth',2,'Color',[0 0.5 1])
hold on
plot(R_vec,N2,'LineWidth',2,'Color',[0 0.5 1],'LineStyle',':')
str = sprintf('Aggr. stats for $z=%g, \\kappa=%g, r=%g, b=%g, \\sigma=%g$ as $R$ varies', z,kappa,r,b,sigma);
title(str,'Interpreter','latex','FontSize',10)
xlabel('$R$','Interpreter','latex')
leg3=legend('$M(R)$','$n(R)$');
set(leg3,'Interpreter','latex');
figure(4)
plot(R_vec,mstar2,'LineWidth',2,'Color',[0 0.6 0],'LineStyle',':')
hold on
plot(R_vec,mbar2,'LineWidth',2,'Color',[0 0.6 0])
str = sprintf('Optimal policy for $z=%g, \\kappa=%g, r=%g, b=%g, \\sigma=%g$ as $R$ varies', z,kappa,r,b,sigma);
title(str,'Interpreter','latex','FontSize',10)
xlabel('$R$','Interpreter','latex')
leg4=legend('$m^{\ast}(R)$','$\bar{m}(R)$');
set(leg4,'Interpreter','latex');
%% 5) Comp statics wrt R (MO)
r=0.005;                                                               % Discount rate of the agent                                                                % Opportunity (flow) cost
b=0.05;                                                                
sigma=0.068;
z=0.17;
kappa=0;
R=0.05;
mstar_init=((3*b*sigma^2)/(4*R))^(1/3);
mbar_init=3*mstar_init;
policy_init=[mstar_init,mbar_init];
R_vec=[0.01:0.005:0.25];
for i=1:length(R_vec)
    R=R_vec(i);
    lambda1=sqrt(2*(r+kappa)/sigma^2);
    lambda2=-sqrt(2*(r+kappa)/sigma^2);
    param_V=[((kappa*b)/r)   R/r    0   R/r    b/r];                       % Parameters given to fun_findV to find V
    param_M=[0                1     0   1       0];                        % Parameters given to fun_findV to find M
    param_n=[kappa            0     0   0       1];     
    fprintf('Calculating statistics and optimal policy for R=%g. \n',R)
    LB=[0,0]; 
    options = optimset('Display','off','TolX',1e-10,'FunValCheck','on','TolCon',1e-5); %Display each iteration
    [policy_opt,fval,EXITFLAG,OUTPUT]=fmincon(@(policy) fun_vstar(policy,param_V),...
    policy_init,[1,-1],0,[],[],LB,[],@constraint,options );          %Constraints: lower bound and mbar>mstar.
    mstar3(i)=policy_opt(1);
    mbar3(i)=policy_opt(2);
    [V3(i),~,~,~,~,~,~,~,~,~,~,~,~]=fun_vstar(policy_opt,param_V);
    [M3(i),~,~,~,~,~,~,~,~,~,~,~,~]=fun_vstar(policy_opt,param_M);
    [N3(i),~,~,~,~,~,~,~,~,~,~,~,~]=fun_vstar(policy_opt,param_n);
end

figure(5)
plot(R_vec,M3,'LineWidth',2,'Color',[0 0.5 1])
hold on
plot(R_vec,N3,'LineWidth',2,'Color',[0 0.5 1],'LineStyle',':')
str = sprintf('Aggr. stats for $z=%g, \\kappa=%g, r=%g, b=%g, \\sigma=%g$ as $R$ varies', z,kappa,r,b,sigma);
title(str,'Interpreter','latex','FontSize',10)
xlabel('$R$','Interpreter','latex')
leg5=legend('$M(R)$','$n(R)$');
set(leg5,'Interpreter','latex');
figure(6)
plot(R_vec,mstar3,'LineWidth',2,'Color',[0 0.6 0],'LineStyle',':')
hold on
plot(R_vec,mbar3,'LineWidth',2,'Color',[0 0.6 0])
str = sprintf('Optimal policy for $z=%g, \\kappa=%g, r=%g, b=%g, \\sigma=%g$ as $R$ varies', z,kappa,r,b,sigma);
title(str,'Interpreter','latex','FontSize',10)
xlabel('$R$','Interpreter','latex')
leg6=legend('$m^{\ast}(R)$','$\bar{m}(R)$');
set(leg6,'Interpreter','latex');

figure(7)
plot(R_vec,M2,'LineWidth',2,'Color',[0 0.5 1])
hold on
plot(R_vec,M3,'LineWidth',2,'Color',[0 0.5 1],'LineStyle',':')
title('Comparative statics of $M$ wrt $R$','Interpreter','latex','FontSize',10)
xlabel('$R$','Interpreter','latex')
ylabel('$M(R)$','Interpreter','latex')
leg7=legend('$\kappa>0$','$\kappa=0$');
set(leg7,'Interpreter','latex');
%% Saving figures
for i=1:1:7
  saveas(figure(i),fullfile('C:\Users\Elia\Dropbox\RoME\MScThesis\Plots\PlotsModel',['compstat' num2str(i) '.png']))
end