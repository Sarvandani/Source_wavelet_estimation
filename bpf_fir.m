function [Dbpf,Data_f_bpf,f]=bpf_fir(Data,dt,N,cut_off)
% This code is written to find the f-x representation of seismic data. The required inputs are:
% Data: complete set of shot gathers
% dt: time sampling interval in seconds
% N: the BPF FIR filter length
% cut_off: a 1x2 vector equal to f_low and f_high cut-off's values
% 
% The outputs are:
% Dbpf: BP filtered Data
% Data_f_bpf: BP filtered Data spectra
% f: a vector containing the frequency values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code if for the book titled: Processing Seismic Reflection
% Data Using MATLAB by Wail A. Mousa & Abdullatif A. Al-Shuhail
% September 2011.

[nt,nx]=size(Data);
nt_fft=2*nt;
Data_f=fx(Data,dt,nt_fft);
B=fir1(N,cut_off*2*dt);
HH=fft(B,nt_fft);
HH=conj(HH)';

Data_f_bpf=zeros(size(Data_f));
for i=1:nx
    Data_f_bpf(:,i)=Data_f(:,i).*HH;
end
Dbpf=real(ifft(Data_f_bpf,[],1));
Dbpf=Dbpf(1+(length(B)-1)/2:nt+(length(B)-1)/2,:);