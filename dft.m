function [ Xk ] = dft( xn , N )
n = [0:1:N-1];
k = [0:1:N-1];
WN = exp(-1i*2*pi/N);
%n*n-tensor construction
nk = n'*k;
WNnk = WN.^nk;
Xk = xn*WNnk;
end

