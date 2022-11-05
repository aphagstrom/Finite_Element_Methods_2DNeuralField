
trimesh(plt1.T,x,y,z)
hold on
plot3(p(1,1),p(2,1),p(3,1),'o','MarkerSize',10,...
    'MarkerEdgeColor','black', 'MarkerFaceColor','black')
axis tight

figure
%R
plot(U(1,1:1000),'red')
title('Mean Firing Rate (R)');
axis tight
xlabel('time steps (dt=0.1)') 
ylabel('R') 

%V
figure
plot(U(np+1,1:1000),'black')
title('Mean Membrane Potential (V)');
axis tight
xlabel('time steps (dt=0.1)') 
ylabel('V') 

%%
W=pi*U(1,1:1000)+1i*U(np+1,1:1000);
Z=(1-conj(W))./(1+conj(W));

figure
plot(1:1000,abs(Z));
title('Within-Population Synchrony (|Z|)');
axis tight
xlabel('time steps (dt=0.1)') 
ylabel('|Z|') 

