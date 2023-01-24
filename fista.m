function [x,J] = fista(y,H,lambda,alpha,Nit)
% [x, J] = fista(y, H, lambda, alpha, Nit)
% L1-regularized signal restoration using the iterated
% soft-thresholding algorithm (FISTA)
% Minimizes J(x) = norm2(y-H*x)^2 + lambda*norm1(x)
% INPUT
% y - observed signal
% H - matrix or operator
% lambda - regularization parameter
% alpha - need alpha >= max(eig(H'*H))
% Nit - number of iterations
% OUTPUT
% x - result of deconvolution
% J - objective function
J = zeros(1, Nit); % Objective function
xkm1 = 0*H'*y; % Initialize xk-1
x = xkm1;
z = x;
T = lambda/(2*alpha);
t = 1.0;
for k = 1:Nit

    Hz = H*z;
    Hxkm1 = H*xkm1;
    
    J(k) = sum(abs(Hxkm1(:)-y(:)).^2) + lambda*sum(abs(xkm1(:)));

    ztmp = z + (H'*(y - Hz))/alpha;
    x = sign(ztmp).*max(0,abs(ztmp)-T);

    a = t-1;
    t = 0.5*(1+sqrt(1+4*t*t));
    a = a / t;
    
    z = x + a*(x-xkm1);
    xkm1 = x;
    
end
