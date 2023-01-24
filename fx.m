function [Data_f,f]=fx(Data,dt,nt_fft)
% This code is written to find the f-x representation of seismic data. The required inputs are:
% Data: complete set of shot gathers
% dt: time sampling interval in seconds
% nt_fft: number of FFT frequency points
% 
% The outputs are:
% Data_f: f-x of Data
% f: a vector containing the frequency values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code if for the book titled: Processing Surface Seismic Reflection
% Data Using MATLAB by Wail A. Mousa & Abdullatif A. Al-Shuhail
% September 2011.

[nt,nx]=size(Data);
if nargin<3
    nt_fft=2*nt;  
end
%f-x data
Data_f=fft(Data,nt_fft,1);
f=linspace(-0.5,0.5,nt_fft)/dt;
