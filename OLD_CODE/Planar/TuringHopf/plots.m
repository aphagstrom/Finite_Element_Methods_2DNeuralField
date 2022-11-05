figure
%R
plot(U(3595,1:9000))
title('Mean Firing Rate (R)');
xlabel('time steps (dt=0.01)') 
ylabel('R') 
axis tight

%V
figure

plot(U(np+3595,1:9000))
title('Mean Membrane Potential (V)') 
ylabel('V') 
axis tight

W=pi*U(3595,1:9000)+1i*U(np+3595,1:9000);
Z=(1-conj(W))./(1+conj(W));

figure
plot(1:9000,abs(Z));
title('Synchrony');
xlabel('time steps (dt=0.01)') 
ylabel('|Z|') 
axis tight

