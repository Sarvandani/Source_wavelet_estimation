function [sx,sy,gx,gy,num_sg,num_trace_per_sg,sz,gz]=extracting_geometry(SegyHeader)
% This code is written to extract the commonly used seismic shot gahter geometries.
% The required inputs is:
% SegyHeader: the seismic data header information
% 
% The outputs are:
% sx: sources x-axis locations
% sy: sources y-axis locations
% gx: receivers x-axis locations
% gy: receivers y-axis locations
% num_sg: is the number of shot gathers
% num_trace_per_sg: number of traces/shot gather
% sz: elevations of the sources
% gz: elevations of the receivers 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code if for the book titled: Processing Surface Seismic Reflection
% Data Using MATLAB by Wail A. Mousa & Abdullatif A. Al-Shuhail
% September 2011.

sx=[SegyHeader.sx];
sy=[SegyHeader.sy];
gx=[SegyHeader.gx];
gy=[SegyHeader.gy];
sz=[SegyHeader.selev];
gz=[SegyHeader.gelev];

sg=[SegyHeader.fldr];
count=1;
l=1;
k=2;
for i=1:length(sg)
    if k<=length(sg) && sg(k-1)==sg(k)
        count=count+1;
    else
        num_trace_per_sg(l)=count;
        num_sg(l)=l;
        count=1;
        l=l+1;
    end
    k=k+1;
end

        
    
