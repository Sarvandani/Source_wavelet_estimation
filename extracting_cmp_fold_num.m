function [cmps,fold_cmp]=extracting_cmp_fold_num(Data,SegyHeader);
% This code is written to extract seismic CMP fold number. The required inputs are:
% Data: complete set of shot gathers
% SegyHeader: Header of the data
% 
% The outputs are:
% cmps: a vector containing the CMP values
% fold_cmp: a vector containing the fold values for each CMP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code if for the book titled: Processing Seismic Reflection
% Data Using MATLAB by Wail A. Mousa & Abdullatif A. Al-Shuhail
% September 2011.

[nt,nx]=size(Data);

for i=1:nx
    XX(i)=SegyHeader(i).cdp;
end

cmps=XX;
cmp_step=cmps(3)-cmps(2);
cmps=min(cmps):cmp_step:max(cmps);
for i=1:length(cmps)
    fold_cmp(i)=length(find(XX==cmps(i)));
end

