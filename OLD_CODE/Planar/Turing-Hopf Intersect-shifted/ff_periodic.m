function [dydt] = ff_periodic(U,param)

np=param.np;
R =  U(1:np);
V =  U(np+1:2*np);
psi = U(2*np+1:3*np);
A1 = U(3*np+1:4*np);
A2 = U(4*np+1:5*np);
A3 = U(5*np+1:6*np);
pp = U(6*np+1:7*np);
g =  U(7*np+1:8*np);

v=param.v;
eta_0=param.eta_0 ;
Delta=param.Delta ;
alfa=param.alfa;
kappaV=param.kappaV;
kappaS=param.kappaS;
tau=param.tau;
M=param.M;
K=param.K;
B=param.B;


% Dxxyy = @(x)-M\(K*x);
Dxxyy = @(u)Dxxyy_periodic(u,K,M,B);
dRdt = (1/tau)*(-kappaV * R + 2 * R .* V + Delta / (tau * pi));
dVdt = (1/tau)*(kappaS*g - (pi * tau * R) .^ 2 + V .^ 2 + eta_0);
d2Rdt2 = (1 / tau) * (-kappaV + 2 * V) .* dRdt + (2 / tau) * R .* dVdt;
A4 = (3/2)*(Dxxyy(A2)) - ((1/v)*dRdt + (1/v^2)*d2Rdt2 - (3/2) * Dxxyy(R));

dydt = [
       dRdt;
    dVdt;
     v * (-psi + A1);
    v * (-A1 + A2 + (3/2)*Dxxyy(psi));
    v * (-A2 + A3);
    v * (-A3 + A4);
     alfa * (-pp + psi);
    alfa * (-g + pp);

];
end