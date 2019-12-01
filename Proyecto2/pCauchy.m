function [pC] = pCauchy( B, g, radio )
% In : B ... (symmetric matrix) approximates the hessian of f in xk
% g ... (vector) gradient of f in xk
% delta ... trust region radius
%
% Out: pC ... The Cauchy point
    pk = g/ norm(g);
    valor = dot(g,B*g);
    if valor > 0
        alphaStar = (norm(g)^3)/(radio*valor);
    else
        alphaStar = 1;
    end
    pC = -radio*min(1,alphaStar)*pk;
    
end

