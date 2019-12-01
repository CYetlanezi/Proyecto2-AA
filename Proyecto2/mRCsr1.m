function [xk, iter] = mRCsr1( f, x0, maxIter, radioMax)
% Trust region method decribed on the previous page
%
% In : f ... (handle) function to be optimized
% x0 ... (vector) initial point
% itmax ... (natural number) upper bound for number of iterations
%
% Out: x ... (vector) last approximation of a stationary point
% iter... (natural number) iterations used
% msg ... (string) message that says whether (or not) a minimum was found


    eta = 0.1;
    r = 1e-6;
    
    radio = radioMax;
    iter = 0;
    xk = x0;
    g = apGrad(f, xk);
    H = speye(length(xk));
    B = H;
    tol = 1e-5; 
    
    while norm(g, 'inf') > tol && iter < maxIter
        %P1
        s = -H*g;
        if dot(s,g) < 0
            if norm(s) > radio
                s = radio*s/norm(s);
            end %else s = s
        else % si no es dir descenso tomamos cauchy
            s = pCauchy( B, g, radio); %checar y probar
        end
         
        %P2
        coc = -( f(xk) - f(xk+s) ) / ( dot(g,s) + 0.5*dot(s,B*s) );
        gnew = apGrad(f, xk + s);
        gamma = gnew - g;
        
        %P3
        if(coc > eta)
            xk = xk + s;
            g = gnew; %
            iter = iter + 1;
            if iter == maxIter
                break
            end
        end % else xk = xk
        
        if coc > 0.75 %P4
            if norm(s) > 0.8*radio
                radio = min(2*radio, radioMax);
            end %else radio = radio
        elseif coc < 0.1 %P5
            radio = 0.5*radio;
            %else radio = radio
        end
            
        v = gamma -B*s;
        %P6
        if abs(dot(v,s)) >= r*norm(v)*norm(s)
            B = B + ( v/dot(v,s) ) * v';
            u = s - H*gamma;
            H = H + ( u/dot(u,gamma) )*u';
        end
    end
end 