function [x,J] = wfista(y,ly,wpow,lambda,alpha,Nit)
%
% L1-regularized signal restoration using the iterated
% soft-thresholding algorithm (FISTA)
% Minimizes J(x) = norm2(y-H*x)^2 + lambda*norm1(x)
%
% adapted to l1 wavelet decon using ffts for matrix applications
%
% INPUT
%
% y - observed signal of length ly
% ensure that ly is good to be used as lfft for ffts (best power of two)
% after correlation with wavelet
%
% wpow(ly) - powerspectrum of wavelet normalized to maximum equal 1
%            real array, twosided for positive as well as negative frequencies
%
% lambda - regularization parameter
% alpha - need alpha >= max(eig(H'*H))
% Nit - number of iterations
%
% OUTPUT
% x - result of L1 deconvolution
% J - objective function
%
J = zeros(1, Nit); % Objective function
xkm1 = zeros(ly,1); % Initialize xk-1
x = zeros(ly,1);
Hz = zeros(ly,1);
ztmp = zeros(ly,1);
z = zeros(ly,1);
Hxkm1 = zeros(ly,1);

x = xkm1;
z = x;
T = lambda/(2*alpha);
t = 1.0;
%

for k = 1:Nit

    %Hz = H*z; % convolve with absorption operator
    
    ftmp = fft(z,ly);
    for n=1:ly
        ftmp(n) = ftmp(n)*wpow(n);
    end
    Hz = real(ifft(ftmp,ly));
    
    %Hxkm1 = H*xkm1; % convolve with absorption operator
    
    ftmp = fft(xkm1,ly);
    for n=1:ly
        ftmp(n) = ftmp(n)*wpow(n);
    end
    Hxkm1 = real(ifft(ftmp,ly));
    
    J(k) = sum(abs(Hxkm1(:)-y(:)).^2) + lambda*sum(abs(xkm1(:)));

    %ztmp = z + (H'*(y - Hz))/alpha; % correlation with absorption operator
    
    for n=1:ly
        ztmp(n) = y(n)-Hz(n);
    end
    ftmp = fft(ztmp,ly);
    for n=1:ly
        ftmp(n) = ftmp(n)*wpow(n);
    end
    ztmp = z + real(ifft(ftmp,ly))/alpha;
    
    x = sign(ztmp).*max(0,abs(ztmp)-T);

    a = t-1;
    t = 0.5*(1+sqrt(1+4*t*t));
    a = a / t;
    
    z = x + a*(x-xkm1);
    xkm1 = x;
    
end
