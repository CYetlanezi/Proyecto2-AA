function [xk, iter] = lineLMBGFS( f, x0, tol, maxiter, m)
    n = length(x0);
    iter = maxiter;
    xk = x0;
    g = apGrad(f, xk);
    H = speye(n);
    S = zeros(n, m);
    Gam = S;
    
    dk = -g;
    
    for k = 1:maxiter
        if norm(g,'inf') <= tol 
            iter = k-1;
            break
        end
%         dk = -H*g;%Paso 1
        assert( dot(dk, g) < 0 )
        [alpha, gnew] = encAlphaLSA( f, xk, dk, g );%paso 2
        
        % Memorizar 
        s   = alpha*dk;
        xk  = xk + s;
        S   = [alpha*dk, S(:, 1:m-1)];
        Gam = [gnew - g, Gam(:, 1:m-1)];
        
        if k < m
            dk = -calcHg(S(:,1:k),Gam(:,1:k),gnew);
        else
            dk = -calcHg(S,Gam,gnew);
        
        end
        g = gnew;
    end
    
end

function[q] = calcHg(S, Gam, gnew)
    % S = {mas nuevo | ...|mas viejo]
    q = gnew;
    m = length(S(1,:));
    irhos = dot(S, Gam, 1);
    alphas = zeros(m,1);
    for i = 1:m
        alphas(i) = dot(S(:,i),	q)/irhos(i);
        q = q - alphas(i)*Gam(:,i);
    end
    
    delta = irhos(1)/dot(Gam(:,i),Gam(:,i));
    q = delta*q;
    
    for i = m:-1:1
        beta = dot(Gam(:,i),q)/irhos(i);
        q = q + (alphas(i) - beta)*S(:,i);
    end
        
end
