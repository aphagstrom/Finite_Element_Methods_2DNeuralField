function [F]  = find_steady_state_NEW(T,param)
R0 = T(1);
V0 = T(2);

%kappaV = param.kappaV;
Delta = param.Delta;
tau = param.tau;
eta_0 = param.eta_0;
maxiter=param.maxiter;
kappaV=param.kappaV;
F=zeros(2,maxiter);

F=[
     -kappaV*R0+2*R0*V0+Delta./(tau*pi);
    
     -(pi * tau * R0) .^ 2 + V0.^ 2 + eta_0];

end

