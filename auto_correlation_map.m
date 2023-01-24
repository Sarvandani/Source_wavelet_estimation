function [Dauto,lags]=auto_correlation_map(Data,max_lag,dt)
% This code is written to obtain the seismic traces auto-correlation map. The required inputs are:
% Data: seismic shot gather(s)
% max_lag: maximum lag value
% dt: time sampling interval
% 
% The output is:
% Ds: the deconvolved seismic data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code if for the book titled: Processing Seismic Reflection
% Data Using MATLAB by Wail A. Mousa & Abdullatif A. Al-Shuhail
% September 2011.

N=round(max_lag/dt);
[nt,nx]=size(Data);
Dauto=zeros(N,nx);
for i=1:nx
    Dauto(:,i)=my_xcorr(Data(:,i)',N)'; 
end
lags=[0:N-1]*dt;
