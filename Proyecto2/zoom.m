function [alphaStar, gnew] = zoom( alpha_lo, alpha_hi, f, xk, dk,gk)
c1 = 1e-4;
c2 = 0.99;
stop = 0;% 0-verdadero 1-falso
maxIter= 200;
k=0;
phi_0 = f(xk);
phi_0_prima = dot( gk, dk);
alphaStar = (alpha_lo + alpha_hi)/2;
gnew = gk;
while stop == 0 && k < maxIter
    alphaj = (alpha_lo + alpha_hi)/2;
    phi_alphaj = f(xk + alphaj*dk);
    phi_alpha_lo = f(xk + alpha_lo*dk);
%     display("while zoom ")
    if (phi_alphaj > (phi_0 + c1*alphaj*phi_0_prima)) ||  ( phi_alphaj >=  phi_alpha_lo )
        alpha_hi = alphaj;
    else
        gnew = apGrad( f, xk + alphaj*dk);
        phi_alphaj_prima = dot( gnew, dk );
        if abs(phi_alphaj_prima) <= -c2*phi_0_prima
            alphaStar = alphaj;
            stop = 1;
            break;
        else
            if phi_alphaj_prima*(alpha_hi - alpha_lo) >= 0
                alpha_hi = alpha_lo;
            end
            alpha_lo = alphaj;
            k = k+1;
        end
    end
end

end