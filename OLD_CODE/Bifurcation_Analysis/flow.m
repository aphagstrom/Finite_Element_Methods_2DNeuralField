clc
clear

tspan = [0 10000]; 
v = 1
eta_0 =3;
Delta = 0.5
alfa = 1
kappaV = 0.5 % Try kappaV=0.5 and then try kappaV=0.6.
kappaS = 1.0
tau = 1
maxiter=200;
axis([0 2.5 0 1])
param.v = v;
param.eta_0 = eta_0;
param.Delta = Delta;
param.alfa = alfa;
param.kappaV = kappaV;
param.kappaS = kappaS;
param.tau = tau;
param.maxiter=maxiter;

F0 = [0.5,0.5]';
options = optimoptions('fsolve','Display','iter'); % Option to display output
y0 = fsolve(@(F)find_steady_state_NEW(F,param),F0,options); % Call solver

ode = @(t,y) neuralmassmodel(t,y,param);
[t,y] = ode45(ode, tspan, y0);

figure
plot(y(:,1),y(:,2),'-o')
axis(2*[-1 1 -1 1])
axis equal
title('Flow at k_v=0.5')
xlabel('R')
ylabel('V')


