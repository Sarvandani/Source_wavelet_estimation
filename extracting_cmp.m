function [D,dt,t,cdp,jj,cmp_offset1]=extracting_cmp(Data,SegyHeader,cmp_num)
% This code is written to extract seismic CMP gathers. The required inputs are:
% Data: complete set of shot gathers
% SegyHeader: Header of the data
% cmp_num: the CMP number
% 
% The outputs are:
% D: The CMP data
% dt: time sampling interval
% t: time vector
% cdp: cdp numbers
% jj: cmp_number locations
% cmp_offset1: cmp offset values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code if for the book titled: Processing Seismic Reflection
% Data Using MATLAB by Wail A. Mousa & Abdullatif A. Al-Shuhail
% September 2011.


[nt,nx]=size(Data);
XX=[SegyHeader.cdp];
cmp_offset=[SegyHeader.offset];

[ii,jj,XX]=find(XX==cmp_num);
cmp_offset1=cmp_offset(jj);
dt = SegyHeader(1).dt/1000/1000;
t=[0:1:nt-1]*dt;
cdp=1:length(jj);
D=Data(:,jj);
