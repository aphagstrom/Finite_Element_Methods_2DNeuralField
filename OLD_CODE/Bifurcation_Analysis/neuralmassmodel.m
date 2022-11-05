function dydt = neuralmassmodel(t,y,param)
v=param.v;
eta_0=param.eta_0;
Delta=param.Delta;
alfa=param.alfa;
kappaV=param.kappaV;
kappaS=param.kappaS;
maxiter=param.maxiter;
tau=param.tau;

dydt = [-kappaV*y(1)+2*y(1)*y(2)+Delta./(tau*pi);
       -(pi * tau * y(1)) .^ 2 + y(2).^ 2 + eta_0];
end