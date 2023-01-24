function [data_out] = inverse_filter(data_in,w_in,w_out,eps)

% data_in - input data
% w_in - input wavelet 
% w_out - desired output wavelet 
% eps - stabilization factor

% number of traces
ntrac = length(data_in(1,:));

% number of time samples
nt = length(data_in(:,1)); 

% number of fft samples
nfft = 2.^nextpow2(nt);

% number of frequencies
nf = nfft/2 + 1;

% fft of the input data
dataf_in = fft(data_in,nfft); 

% keep only one side of the fft
dataf_in = dataf_in(1:nf,:);

% fft of the input wavelet
wf_in = fft(w_in,nfft); wf_in = wf_in(1:nf);

% fft of the output wavelet
wf_out = fft(w_out,nfft); wf_out = wf_out(1:nf);

% cross-correlate the input data with the input wavelet
dataf_ac = dataf_in.*conj(repmat(wf_in,1,ntrac));

% auto-correlate the input wavelet
denom = (abs(wf_in.^2));

% water level for stabilization
wl = eps.^2.*max(denom); 
denom(denom<wl) = wl;

% use repmat to make denom and wf_out the right size for matrix
% multiplication/division
denom = repmat(denom,1,ntrac);
wf_out = repmat(wf_out,1,ntrac);

% divide by the auto-correlation and multiply by the output wavelet
dataff_out = wf_out.*dataf_ac./denom;

% put the two sides of the fft together
dataf_out = zeros(nfft,ntrac);
dataf_out(1:nf,:) = dataff_out;
dataf_out(nf+1:end,:) = conj(dataff_out(nf-1:-1:2,:));

% inverse fourier transform
data_out = ifft(dataf_out);

% keep the first nt time samples to match the input data
data_out = data_out(1:nt,:);
