function[out] = ifft_dfh_2 (in,nf,n,nt)


% inputs:
% in – data in the f-x domain (N-by-M, where N is the number of frequencies from 1 to nyquist, and M is the number of traces); 
% nf  - number of frequency samples from 1 to nyquist 
% n – total number of frequency samples
% nt – desired number of output time samples
% 
% tI always use the following realtionship between nf, n, and nt
% n = 2.^nextpow2(nt); 
% nf = n/2 + 1;
% corresponding frequency axis should be:
%  freq = [0:nf-1].*df; 
% with:
% df = 1/(n*dt);


nx = length(in(1,:));

yf = zeros(n,nx);
out = zeros(nt,nx);
% for ii = 1:nx
%     clear yf
%     yf(1:nx,1:n) = 0;
    yf(2:nf,1:nx) = in(2:nf,1:nx);
    
%     for ifreq = nf+1:n
        yf(nf+1:n,1:nx) = conj(yf((nf-1):-1:2,1:nx));
%     end
    
% end

    
    yt = real(ifft(yf));
%     out(nt:end,ii) = yt(1,1:nt);
    out(1:nt,:) = yt(1:nt,:);
end