function [alphaStar, gnew] = encAlphaLSA( f, xk, dk, gk )
c1 = 1e-4;
c2 = 0.99;
alphai_ant = 0; %alpha0
alphai = 1; %alpha1
alphamax = 128;
i = 1;
maxIter = 200;
stop = 0;% 0-verdadero 1-falso
phi_0 = f(xk);
phi_0_prima = dot(gk,dk);
k =0;
while stop == 0 && k <  maxIter
    phi_alphai = f(xk + alphai*dk); 
    phi_alphai_ant = f(xk + alphai_ant*dk); 
    
    if (phi_alphai > phi_0 + c1*alphai*phi_0_prima) ||  ( phi_alphai >=  phi_alphai_ant && i > 1)
        [alphaStar, gnew] = zoom(alphai_ant,alphai,f, xk, dk,gk); 
%         alphaStar = zoom(alphai_ant,alphai,f, xk, dk,gk);
%         gnew = apGrad(f,xk + alphaStar*dk);
        stop = 1;
        break;
    else
        gnew = apGrad(f, xk + alphai*dk);
        phi_alphai_prima = dot( gnew, dk);
        if abs(phi_alphai_prima) <= -c2*phi_0_prima
            alphaStar = alphai;
            stop = 1;
            break;
        else
            if phi_alphai_prima >= 0
               [alphaStar, gnew] = zoom(alphai_ant,alphai,f, xk, dk,gk); 
%                 alphaStar = zoom(alphai,alphai_ant,f, xk, dk,gk);
%                 gnew = apGrad(f,xk + alphaStar*dk);
                stop = 1;
                break;
            else
%                 alphai_ant = alphai;
                alphai = min(alphamax, 2*alphai);
                i = i+1;
            end
        end
    end
end
end