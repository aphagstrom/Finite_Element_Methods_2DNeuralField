
trimesh(plt1.T,x,y,z)
hold on
plot3(p(1,1),p(2,1),p(3,1),'o','MarkerSize',10,...
    'MarkerEdgeColor','black', 'MarkerFaceColor','black')
axis tight

figure
%R
plot(U(1,1:18000),'red')
title('Mean Firing Rate (R) at k_v=0.6');
axis tight
xlabel('time steps (dt=0.1)') 
ylabel('R') 

%V
figure

plot(U(np+1,1:18000),'black')
title('Mean Membrane Potential (V) at k_v=0.6');
axis tight
xlabel('time steps (dt=0.1)') 
ylabel('V') 

W=pi*U(1,1:18000)+1i*U(np+1,1:18000);
Z=(1-conj(W))./(1+conj(W));

figure
plot(1:18000,abs(Z));
title('Within-Population Synchrony (|Z|) and k_v=0.6');
axis tight
xlabel('time steps (dt=0.1)') 
ylabel('|Z|') 

