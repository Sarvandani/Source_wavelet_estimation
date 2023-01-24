function Ds=spiking_decon(Data,max_lag,mu,dt)
% This code is written to perform spiking deconvolution for seismic traces.
% It uses the auto-correlation generated based on the whole record time.
%
% The required inputs are:
% Data: seismic shot gather(s)
% max_lag: maximum lag value(Note that the min lag is set to be the first sample)
% mu: White noise percentage
% dt: time sampling interval
% 
% The output is:
% Ds: the deconvolved seismic data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code if for the book titled: Processing Surface Seismic Reflection
% Data Using MATLAB by Wail A. Mousa & Abdullatif A. Al-Shuhail
% September 2011.

N=round(max_lag/dt);
p_noise=mu/100;
[nt,nx]=size(Data);
Dauto=zeros(N,nx);
Ds=zeros(nt+N-1,nx);
Rxd=zeros(1,N);

for i=1:nx
    Dauto(:,i)=my_xcorr(Data(:,i)',N)'; 
    Rxd=Rxd+my_xcorr(Data(:,i)',N,1);
end
Rxd=abs(Rxd);

DDauto=sum(Dauto,2);
Rxx=toeplitz(DDauto);
Rxx(1,1)=Rxx(1,1)*(p_noise);


h_opt=(inv(Rxx+eps))*Rxd';
for i=1:nx
    Ds(:,i)=conv(Data(:,i),h_opt);
end

Ds=Ds(1:nt,:);