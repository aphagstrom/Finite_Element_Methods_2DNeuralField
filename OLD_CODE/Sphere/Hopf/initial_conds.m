function[k1,k2,init] =initial_conds(nterms,kc,phi_vec,c_vec,x,y)
 init = zeros(1,length(x));
for ii=1:nterms
        k1= kc*cos(phi_vec(ii));
        k2=  kc* sin(phi_vec(ii));
        init=init+c_vec(ii).*exp(1i.*k1.*x + 1i.*k2.*y) + conj(c_vec(ii).*exp(-1i.*k1.*x-1i.*k2.*y));
end

end