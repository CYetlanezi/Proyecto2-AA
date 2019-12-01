function [f, Df, Hf] = fDixmaanaG()
  
    f = @(x)  DixmaanaG(x);
    
    Df = @(x)  apGrad(f,x);
    
    Hf = @(x)  apHess(f,x);
end
