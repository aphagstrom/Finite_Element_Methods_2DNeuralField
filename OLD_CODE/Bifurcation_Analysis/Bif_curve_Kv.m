clc
clear

F0 = [0.5,0.5]';
eta_0 =3
Delta = 0.5
alfa = 1
kappaV = 0
kappaS = 1
tau = 1
v=1


maxiter=200;
param.eta_0 = eta_0;
param.Delta = Delta;
param.alfa = alfa;
param.kappaV = kappaV;
param.kappaS = kappaS;
param.tau = tau;
parameter.v=v;
param.maxiter=maxiter;
options = optimoptions('fsolve','Display','iter'); % Option to display output
C=zeros(3,maxiter);
eigenA=zeros(2,maxiter);
AA=zeros(2,maxiter);
expr=zeros(1,maxiter);
for ITER= 1 : maxiter   
  F(1:2,ITER) = fsolve(@(F) find_steady_state_NEW(F, param), F0, options);
  param.kappaV=0.01+param.kappaV;
  F0=F(1:2,ITER);
  [F(1:2,ITER),fval,exitflag,output,jacobian]  = fsolve(@(F)find_steady_state_NEW(F,param),F0,options);
  eigenA(1:2,ITER)=eig(jacobian);
  if real(eigenA(1,ITER))>0
  scatter(param.kappaV,F(1,ITER),'black')
  elseif real(eigenA(1,ITER))<0
  scatter(param.kappaV,F(1,ITER),'red')    
  end
  %axis([0 2.5 0 1])
  axis tight
  xlabel('k_v')
  ylabel('R')
  title('k_v versus steady state')
  hold on
end




 
