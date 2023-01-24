function [D,dt,t,offset]=extracting_shots_cmp(Data,SegyHeader,shot_num,p)
% This code is written to extract a seismic shot gahter or a group of shot
% gathers. The required inputs are:
% Data: complete set of shot gathers
% SegyHeader: the seismic data header information
% shot_num: the shot gather number or a group of shot gathers
% p: is 0 for offset (only for a sigle shot gahter) or 1 for trace numbers
% 
% The outputs are:
% D: the required shot gahter(s)
% dt: time sampling interval in seconds
% dx: spatial sampling interval in meters or feet depending on the data
% offset: is the data x-axis vector containing (depending on the p value above) either:
% (a) the offset values (only if the user requests a single shot gather)
% (b) the traces' numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code if for the book titled: Processing Surface Seismic Reflection
% Data Using MATLAB by Wail A. Mousa & Abdullatif A. Al-Shuhail
% September 2011.

[nt,nx]=size(Data);
XX=[SegyHeader.fldr];
dt = SegyHeader(1).dt/1000/1000;
t=[0:1:nt-1]*dt;

    [ii,jj]=find(XX==shot_num);
    ll=1;
    for i=min(jj):max(jj)
        offset(ll)=SegyHeader(i).offset;
        ll=ll+1;
    end
    
    D=Data(:,jj);
    [nt_D,nx_D]=size(D);

