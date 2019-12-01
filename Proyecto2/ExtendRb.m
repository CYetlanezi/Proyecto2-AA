function y = ExtendRb(x)
n = length(x);
n = n/2;
suma = 0;
for j = 1:n
    suma = suma + 100* (x(2*j) - (x( 2*j - 1))^2 )^2 + (1 - x(2*j -1) )^2;
end
y = suma;
