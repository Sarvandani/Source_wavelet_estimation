function [Data_fk,f,kx]=fk(Data,dt,dx,nt_fft,nx_fft)
% This code is written to find the f-x representation of seismic data. The required inputs are:
% Data: complete set of shot gathers
% dt: time sampling interval in seconds
% nt_fft: number of FFT frequency points
% nx_fft: number of FFT wavenumber points
%
% The outputs are:
% Data_fk: f-kx of Data
% f: a vector containing the frequency values
% kx: a vector containing the wavenumber values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code if for the book titled: Processing Surface Seismic Reflection
% Data Using MATLAB by Wail A. Mousa & Abdullatif A. Al-Shuhail
% September 2011.

[nt,nx]=size(Data);
if nargin<4
    nt_fft=2*nt;
    nx_fft=2*nx;
end
%f-x data
Data_f=fft(Data,nt_fft,1);
Data_fk=fft(Data_f,nx_fft,2);
f=linspace(-0.5,0.5,nt_fft)/dt;
kx=linspace(-0.5,0.5,nx_fft)/dx;
