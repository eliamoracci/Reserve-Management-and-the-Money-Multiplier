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
z=0.22;
mstar_init=((3*b*sigma^2)/(4*R))^(1/3);
mbar_init=3*mstar_init;
policy_init=[mstar_init,mbar_init];
kappa_vec=[0:0.01:0.16 0.17:0.001:0.18 0.19 0.2];
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
    
    n=3;
mstar_init_vec=mstar_init*0.8:(mstar_init*1.2-mstar_init*0.8)/(n-1):mstar_init*1.2;
mbar_init_vec=mbar_init*0.8:(mbar_init*1.2-mbar_init*0.8)/(n-1):mbar_init*1.2;
policy_init_vec=cell(n,n);
for ii=1:n
    for j=1:n
        policy_init_vec(j,ii)={[mstar_init_vec(j) mbar_init_vec(ii)]};
    end
end
fval_vec=cell(n,n);
policy_opt_vec=cell(n,n);
options = optimset('Display','off'); %Display each iteration or not
for j=1:n
    for ii=1:n
        [policy_opt_a,fval,EXITFLAG,OUTPUT]=fmincon(@(policy) fun_vstar(policy,param_V),...
    cell2mat(policy_init_vec(j,ii)),[1,-1],0,[],[],LB,[],@constraint,options );    
    policy_opt_vec(j,ii)={[policy_opt_a(1) policy_opt_a(2)]};
    fval_vec(j,ii)={fval};
    end
end
fval_vec=cell2mat(fval_vec);
[x,y]=find(fval_vec==min(min(fval_vec)));
policy_opt_final=cell2mat(policy_opt_vec(x,y));
mstar_opt_final=policy_opt_final(1);
mbar_opt_final=policy_opt_final(2);
tol=0.0001;
if abs(mstar_opt_final-mstar1(i))<tol && abs(mbar_opt_final-mbar1(i))<tol
    disp('I found the optimal policy!')
else
    mstar1(i)=mstar_opt_final;
    mbar1(i)=mbar_opt_final;
    policy_opt=[mstar1(i) mbar1(i)];
    fprintf('I had found a local minimum, now I fixed it: mstar=%g, mbar=%g. \n',mstar1(i),mbar1(i));
end
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
%% Comp statics wrt R 
r=0.005;                                                               % Discount rate of the agent                                                                % Opportunity (flow) cost
b=0.05;                                                                
sigma=0.05;
z=0.22;
kappa=0.1;
R=0.05;
mstar_init=((3*b*sigma^2)/(4*R))^(1/3);
mbar_init=3*mstar_init;
policy_init=[mstar_init,mbar_init];
R_vec=[0.01:0.01:0.2];
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
        
    n=3;
mstar_init_vec=mstar_init*0.8:(mstar_init*1.2-mstar_init*0.8)/(n-1):mstar_init*1.2;
mbar_init_vec=mbar_init*0.8:(mbar_init*1.2-mbar_init*0.8)/(n-1):mbar_init*1.2;
policy_init_vec=cell(n,n);
for ii=1:n
    for j=1:n
        policy_init_vec(j,ii)={[mstar_init_vec(j) mbar_init_vec(ii)]};
    end
end
fval_vec=cell(n,n);
policy_opt_vec=cell(n,n);
options = optimset('Display','off'); %Display each iteration or not
for j=1:n
    for ii=1:n
        [policy_opt_a,fval,EXITFLAG,OUTPUT]=fmincon(@(policy) fun_vstar(policy,param_V),...
    cell2mat(policy_init_vec(j,ii)),[1,-1],0,[],[],LB,[],@constraint,options );    
    policy_opt_vec(j,ii)={[policy_opt_a(1) policy_opt_a(2)]};
    fval_vec(j,ii)={fval};
    end
end
fval_vec=cell2mat(fval_vec);
[x,y]=find(fval_vec==min(min(fval_vec)));
policy_opt_final=cell2mat(policy_opt_vec(x,y));
mstar_opt_final=policy_opt_final(1);
mbar_opt_final=policy_opt_final(2);
tol=0.0001;
if abs(mstar_opt_final-mstar2(i))<tol && abs(mbar_opt_final-mbar2(i))<tol
    disp('I found the optimal policy!')
else
    mstar2(i)=mstar_opt_final;
    mbar2(i)=mbar_opt_final;
    policy_opt=[mstar2(i) mbar2(i)];
    fprintf('I had found a local minimum, now I fixed it: mstar=%g, mbar=%g. \n',mstar1(i),mbar1(i));
end
    
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

%% Saving figures
for i=1:1:4
  saveas(figure(i),fullfile('C:\Users\Elia\Dropbox\RoME\MScThesis\Plots\PlotsModel',['compstat' num2str(i) '.png']))
end