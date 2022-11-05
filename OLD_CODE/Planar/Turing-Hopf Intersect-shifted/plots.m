trimesh(DT,x,y)
hold on
plot(p(1,1677),p(2,1677),'o','MarkerSize',10,...
    'MarkerEdgeColor','black', 'MarkerFaceColor','black')
axis tight

figure
%R
plot(U(1677,1:3994),'red')
title('Mean Firing Rate (R)');
axis tight
xlabel('time steps (dt=0.1)') 
ylabel('R') 
axis tight
%V
figure

plot(U(np+1677,1:3994),'black')
title('Mean Membrane Potential (V)');
axis tight
xlabel('time steps (dt=0.1)') 
ylabel('V') 
axis tight

W=pi*U(1677,1:3994)+1i*U(np+1677,1:3994);
Z=(1-conj(W))./(1+conj(W));

figure
plot(1:3994,abs(Z));
title('Within-Population Synchrony |Z|');
axis tight
xlabel('time steps (dt=0.1)') 
ylabel('|Z|') 


