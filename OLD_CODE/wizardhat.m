x=linspace(-10,10,500);
y = 1/(2*pi).*(abs(x)./2-1).*exp(-abs(x));
p=plot(x,y,'r');

%p(1).Marker='*';
title('Wizard Hat Function w(r)')
xlabel('r') 
ylabel('w') 
p.LineWidth=4