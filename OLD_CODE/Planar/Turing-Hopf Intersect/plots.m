figure

%R
plot(U(3595,1:1800))
%V
figure
plot(U(np+3595,1:1800))

W=pi*U(3595,1:1800)+1i*U(np+3595,1:1800);
Z=(1-conj(W))./(1+conj(W));

figure
plot(1:1800,abs(Z));
