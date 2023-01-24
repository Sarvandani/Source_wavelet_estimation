function [Dstacked,t,cmp_num]=sstack(Dsort,Hsort)
% This code is written to perform stacking for seismic CMP gathers.
% The required inputs are:
% Dsort: the sorted data
% Hsort: the sorted header
% 
% The outputs are:
% Dstacked: the stacked seismic data
% t: time axis values
% cmp_num: a vector containing the corresponding cmp values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code if for the book titled: Processing Surface Seismic Reflection
% Data Using MATLAB by Wail A. Mousa & Abdullatif A. Al-Shuhail
% September 2011.

cmp_num=extracting_cmp_fold_num(Dsort,Hsort);
%Stacking
for kkk=1:length(cmp_num)
    disp(['I am CMP number: ',num2str(cmp_num(kkk)),''])
    [Da,dt,t]=extracting_cmp(Dsort,Hsort,cmp_num(kkk));
    Dstacked(:,kkk)=sum(Da,2);
       
end