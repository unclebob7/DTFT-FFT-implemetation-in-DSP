function [ y ] = cirshift( x , m , N )
if length(x)>N
    error('N should be larger or equal to x(n)-length');
end
x = [x , zeros(1,N-length(x))];
n = 0:(N-1);
n = mod(n-m , N);
y = x(n+1);
end

