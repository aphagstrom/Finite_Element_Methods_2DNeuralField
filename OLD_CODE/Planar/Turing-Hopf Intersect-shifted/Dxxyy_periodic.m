function v = Dxxyy_periodic(u,K,M,B)

vl = -[M B; B' zeros(size(B,2))] \ [K*u; zeros(size(B,2),1)];
v = vl(1:length(u));