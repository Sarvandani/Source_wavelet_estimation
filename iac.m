function D=iac(Data,t,pow,T)
% This code is written to perform indepedent amplitude correction for seismic traces. The required inputs are:
% Data: seismic shot gather(s)
% t: the time axis values
% pow: the power
% T: is 0 for the power of time correction or 1 for exponential gain
% function correction
% 
% The output is:
% D: the amplitude corrected data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code if for the book titled: Processing Surface Seismic Reflection
% Data Using MATLAB by Wail A. Mousa & Abdullatif A. Al-Shuhail
% September 2011.

[nt,nx]=size(Data);
D=zeros(nt,nx);
if T==0
    t=(t').^pow;
    for i=1:nx
        D(:,i)=Data(:,i).*t;
    end

elseif T==1
    t=exp((t').^pow);
    for i=1:nx
        D(:,i)=Data(:,i).*t;
    end
end
    
    