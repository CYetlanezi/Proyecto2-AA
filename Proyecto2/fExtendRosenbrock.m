function [f, Df, Hf] = fExtendRosenbrock()
  
    f = @(x)  ExtendRb(x);
    
    Df = @(x)  apGrad(f,x);
    
    Hf = @(x)  apHess(f,x);
end
