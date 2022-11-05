
alpha=0.25;
x=linspace(-0.5,50);
y=alpha^2.*x.*exp(-alpha.*x).*heaviside(x);
p=plot(x,y)
p.LineWidth=4;
title('Synaptic Function s(t)')
xlabel('t') 
ylabel('s') 

