function y = DixmaanaG(x)
n = length(x);
m = floor(n/3);
alpha = 1;
beta = 0.125;
gamma = 0.125;
delta = 0.125;
k1 = 1;
k2 = 0;
k3 = 0;
k4 = 1;

suma1 = 0;
suma2 = 0;
suma3 = 0;
suma4 = 0;

for j = 1:n
    suma1 = suma1 + alpha *((x(j))^2)*(j/n)^k1;
    if j < n
        suma2 = suma2 + beta *((x(j))^2)*( x(j+1) + (x(j+1))^2 )^2*(j/n)^k2 ;
    end
    if j <= 2*m 
        suma3 = suma3 + gamma *( ( x(j) )^2 )*( (x(j+m))^4 )*(j/n)^k3;
    end
    if j <= m
        suma4 = suma4 + delta *( x(j) )*( x( j + 2*m ))*(j/n)^k4;
    end
end

y = 1 + suma1 + suma2 + suma3 + suma4;
