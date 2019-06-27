%This code solves the augmented Miller-Orr model where the law of motion
%for the stock of the liquid asset is dm(t)=\sigma dW(t) -zdN(t), where
%W(t) is a Wiener process and N(t) is a Poisson process with rate of
%arrival \kappa and negative jump size z.

%In particular, the code:
% 1) Finds the optimal policy [mstar,mbar] for a given initial condition
%    and tries other IC to see if it's a global or local minimum.
% 2) Plots the value function V(m), the n. adj. N(m) and the money holdings
%    M(m) and derives average N and M. It also plots the empirical
%    distribution h(m) deriving from simulating the process subject to the
%    optimal control policy.
% 3) It plots the minimum of the value function V(mstar,mstar,mbar) for different
%    combinations of policy parameters [mstar,mbar], so that one is able to
%    visualize in three dimensions if there are several local minima, if
%    the global minimum is where part 1) indicated, and so on. For the
%    moment, it is plotted for every mstar and mbar in [0,1] because for a
%    large set of parameters the minimum is in this interval.
% 4) It computes again the optimal policy starting from n different initial
%    guesses; tipically different initial guesses will result into
%    different solutions found (multiple local minima). Of the solutions
%    found, I check that the best (the one with the lower V(mstar)) is the
%    one I found at part (1) (starting from the classical Miller-Orr
%    solution as a guess). If it's not the case, I must dismiss what I
%    found in part (2).
%
%This code exploits several functions:
% 1) fun_Vstar (which exploits fsolve and fun_findV)
% 2) fun_findV (which exploits solvesys and solvesys0)
% 3) simulate  (needed to simulate a controlled jump diffusion)
% 4) constraint (needed to impose the smooth pasting condition when using
%                fmincon in part (1) and (4).)
% 5) fun_F (needed to get functional forms for V(m),N(m) and M(m))
%% Housekeeping
clear all
close all
%% Defining variables
global r R b kappa z sigma lambda1 lambda2 param_V param_M param_n
r=0.005;                                                               % Discount rate of the agent
R=0.05;                                                                % Opportunity (flow) cost
b=0.05;                                                                
z=0.05;
kappa=0.0001;
sigma=0.05;
lambda1=sqrt(2*(r+kappa)/sigma^2);
lambda2=-sqrt(2*(r+kappa)/sigma^2);
param_V=[((kappa*b)/r)   R/r    0   R/r    b/r];                       % Parameters given to fun_findV to find V
param_M=[0                1     0   1       0];                        % Parameters given to fun_findV to find M
param_n=[kappa            0     0   0       1];                        % Parameters given to fun_findV to find N
%% 1) Finding the optimal policy [mstar,mbar]

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
disp('Checking global optimality of the solution...')
n=3;
mstar_init_vec=mstar_init*0.8:(mstar_init*1.2-mstar_init*0.8)/(n-1):mstar_init*1.2;
mbar_init_vec=mbar_init*0.8:(mbar_init*1.2-mbar_init*0.8)/(n-1):mbar_init*1.2;
policy_init_vec=cell(n,n);
for i=1:n
    for j=1:n
        policy_init_vec(j,i)={[mstar_init_vec(j) mbar_init_vec(i)]};
    end
end
fval_vec=cell(n,n);
policy_opt_vec=cell(n,n);
options = optimset('Display','off'); %Display each iteration or not
for j=1:n
    for i=1:n
        [policy_opt_a,fval,EXITFLAG,OUTPUT]=fmincon(@(policy) fun_vstar(policy,param_V),...
    cell2mat(policy_init_vec(j,i)),[1,-1],0,[],[],LB,[],@constraint,options );    
    policy_opt_vec(j,i)={[policy_opt_a(1) policy_opt_a(2)]};
    fval_vec(j,i)={fval};
    end
end
fval_vec=cell2mat(fval_vec);
[x,y]=find(fval_vec==min(min(fval_vec)));
policy_opt_final=cell2mat(policy_opt_vec(x,y));
mstar_opt_final=policy_opt_final(1);
mbar_opt_final=policy_opt_final(2);
tol=0.0001;
if abs(mstar_opt_final-mstar_opt)<tol && abs(mbar_opt_final-mbar_opt)<tol
    disp('I found the optimal policy!')
else
    mstar_opt=mstar_opt_final;
    mbar_opt=mbar_opt_final;
    policy_opt=[mstar_opt mbar_opt];
    fprintf('I had found a local minimum, now I fixed it: mstar=%g, mbar=%g. \n',mstar_opt,mbar_opt);
end

%% 2) Plotting value function, avg. n. adj., avg. money h. and empirical h(m)
disp('Making plots of V(m), M(m), N(m) and h(m).')
% Plotting the value function found
n_eval=1000;
m_vec=0:mbar_opt/(n_eval-1):mbar_opt;
[Vstar,B001,B002,~,~,A0,D0,A,D,B01,B02,B1,B2]=fun_vstar(policy_opt,param_V);
[valuef]=fun_F(B001,B002,A0,D0,A,D,B01,B02,B1,B2,m_vec);
value_function=valuef;
figure(1)
plot(m_vec,valuef,'LineWidth',2);
axis([min(m_vec)-0.1*max(m_vec) 1.1*max(m_vec) ...
    min(valuef)-(0.1)*(max(valuef)-min(valuef)) ...
    max(valuef)+(0.1)*(max(valuef)-min(valuef)) ]);
grid on
hline(Vstar,'--r')
vline(0,'--r');
vline(mstar_opt,'r'); 
vline(mbar_opt,'--r'); 
t0=text(mstar_opt+(0.01)*(max(valuef)-min(valuef)),...
    min(valuef)-(0.02)*(max(valuef)-min(valuef)),...
    '$m^{\ast}$','Interpreter','latex');
set(t0,'Fontsize',11)
t1=text(mbar_opt+(0.03)*(max(valuef)-min(valuef)),...
    min(valuef)-(0.06)*(max(valuef)-min(valuef)),...
    '$\bar{m}$','Interpreter','latex');
set(t1,'Fontsize',11)
xlabel('$m$','Interpreter','latex')
ylabel('$V(m,m^{\ast},\bar{m})$','Interpreter','latex')
str = sprintf('Value function for $\\kappa=%g, z=%g, r=%g, R=%g, b=%g, \\sigma=%g$', kappa,z,r,R,b,sigma);
title(str,'Interpreter','latex')

figure(6)
plot(m_vec,valuef,'LineWidth',2);
axis([0.05 0.2 ...
    min(valuef)-(0.03)*(max(valuef)-min(valuef)) ...
    min(valuef)+(0.1)*(max(valuef)-min(valuef)) ]);
grid on
hline(Vstar,'--r')
vline(0,'--r');
vline(mstar_opt,'r'); 
vline(mbar_opt,'--r'); 
t0=text(mstar_opt+(0.01)*(max(valuef)-min(valuef)),...
    min(valuef)-(0.02)*(max(valuef)-min(valuef)),...
    '$m^{\ast}$','Interpreter','latex');
set(t0,'Fontsize',11)
t1=text(mbar_opt+(0.03)*(max(valuef)-min(valuef)),...
    min(valuef)-(0.06)*(max(valuef)-min(valuef)),...
    '$\bar{m}$','Interpreter','latex');
set(t1,'Fontsize',11)
xlabel('$m$','Interpreter','latex')
ylabel('$V(m,m^{\ast},\bar{m})$','Interpreter','latex')
title('Detail','Interpreter','latex')

%% Plot the average money holdings
r=r/10000;                                                             %I take the 'limit' of r to 0 to compute average statistics.
[M,B001,B002,~,~,A0,D0,A,D,B01,B02,B1,B2]=...                          %I exploit param_M to change the ODEs system and find M.
    fun_vstar(policy_opt,param_M);
[valuef]=fun_F(B001,B002,A0,D0,A,D,B01,B02,B1,B2,m_vec);
avg_money_function=valuef;
figure(2)
plot(m_vec,valuef,'LineWidth',2);
axis([min(m_vec)-0.1*max(m_vec) 1.1*max(m_vec) ...
    min(valuef)-(0.1)*(max(valuef)-min(valuef)) ...
    max(valuef)+(0.1)*(max(valuef)-min(valuef)) ]);
grid on
hline(M,'--r')
vline(0,'--r');
vline(mstar_opt,'r'); 
vline(mbar_opt,'--r');  
t0=text(mstar_opt+(0.03)*(max(valuef)-min(valuef)),...
    min(valuef)-(0.06)*(max(valuef)-min(valuef)),...
    '$m^{\ast}$','Interpreter','latex');
set(t0,'Fontsize',11)
t1=text(mbar_opt+(0.03)*(max(valuef)-min(valuef)),...
    min(valuef)-(0.06)*(max(valuef)-min(valuef)),...
    '$\bar{m}$','Interpreter','latex');
set(t1,'Fontsize',11)
xlabel('$m$','Interpreter','latex')
ylabel('$M(m,m^{\ast},\bar{m})$','Interpreter','latex')
str = sprintf('Average exc. res. for $\\kappa=%g, z=%g, r=%g, R=%g, b=%g, \\sigma=%g$',...
    kappa,z,r,R,b,sigma);
title(str,'Interpreter','latex')

% Plot the average number of adjustments
[n_adj,B001,B002,~,~,A0,D0,A,D,B01,B02,B1,B2]=...                      %I exploit param_n to change the ODEs system and find N.
    fun_vstar(policy_opt,param_n);  
[valuef]=fun_F(B001,B002,A0,D0,A,D,B01,B02,B1,B2,m_vec);
avg_adj_function=valuef;
figure(3)
plot(m_vec,valuef,'LineWidth',2);
axis([min(m_vec)-0.1*max(m_vec) 1.1*max(m_vec) ...
    min(valuef)-(0.1)*(max(valuef)-min(valuef)) ...
    max(valuef)+(0.1)*(max(valuef)-min(valuef)) ]);
grid on
hline(n_adj,'--r')
vline(0,'--r');
vline(mstar_opt,'r'); 
vline(mbar_opt,'--r'); 
t0=text(mstar_opt+(0.03)*(max(valuef)-min(valuef)),...
    min(valuef)-(0.06)*(max(valuef)-min(valuef)),...
    '$m^{\ast}$','Interpreter','latex');
set(t0,'Fontsize',11)
t1=text(mbar_opt+(0.03)*(max(valuef)-min(valuef)),...
    min(valuef)-(0.06)*(max(valuef)-min(valuef)),...
    '$\bar{m}$','Interpreter','latex');
set(t1,'Fontsize',11)
xlabel('$m$','Interpreter','latex')
ylabel('$n(m,m^{\ast},\bar{m})$','Interpreter','latex')
str = sprintf('Average n. adj. for $\\kappa=%g, z=%g, r=%g, R=%g, b=%g, \\sigma=%g$',...
    kappa,z,r,R,b,sigma);
title(str,'Interpreter','latex')

r=r*10000;
path=simulate(policy_opt,50000);               
figure(4)
histogram(path,100,'Normalization','pdf','DisplayStyle','bar','FaceColor',[0.3 0.3 0.3],'FaceAlpha',0.2)
str = sprintf('Simulated PDF of $m$ for $\\kappa=%g, z=%g, r=%g, R=%g, b=%g, \\sigma=%g$', kappa,z,r,R,b,sigma);
title(str,'Interpreter','latex','Fontsize',9)
vline(mstar_opt,'r'); 
vline(mbar_opt,'--r');
xlabel('$m$','Interpreter','latex')
ylabel('$$\hat{h}(m)$$','Interpreter','latex')
xticks([0 mstar_opt mbar_opt])
xticklabels({'0','$m^{\ast}$','$\bar{m}$'})
%% 3) Plotting V^{m^{\ast},m^{\ast},\bar{m}}
disp('Plotting V^{m^{\ast},m^{\ast},\bar{m}}.')
n=40;
mstar_vec=0:1/(n-1):1;
mbar_vec=0:1/(n-1):1;
vstar_mat=zeros(n,n);
for i=1:length(mstar_vec)
    for j=1:length(mbar_vec)
        policy=[mstar_vec(i) mbar_vec(j)];
        if mbar_vec(j)>mstar_vec(i)
            [vstar_mat(i,j),~,~,~,~,~,~,~,~,~,~,~,~]=fun_vstar(policy,param_V);
        else
            vstar_mat(i,j)=NaN;
        end
    end
end
figure(5)
surf(mbar_vec,mstar_vec,vstar_mat);
view(225,50)
camlight
lighting gouraud
xlabel('$\bar{m}$','Interpreter','latex')
ylabel('$m^{\ast}$','Interpreter','latex')
hline(mstar_opt,'b--')
vline(mbar_opt,'b--')
zlabel('$V(m^{\ast},m^{\ast},\bar{m})$','Interpreter','latex')
str = sprintf('Value function for $z=%g, r=%g, R=%g, b=%g, \\sigma=%g, \\kappa=%g$', z,r,R,b,sigma,kappa);
title(str,'Interpreter','latex','FontSize',10)
%% Saving figures
for i=1:1:6
  saveas(figure(i),fullfile('C:\Users\Elia\Dropbox\RoME\MScThesis\Plots\PlotsModel',['figure' num2str(i) '.png']))
end

