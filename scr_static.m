function Dsort=scr_static(Dsort,Hsort,cmp_start,cmp_end,maxlags)
% This code is written to perform static correction on the NMO CMP gathers.
% The required inputs are:
% Dsort: the sorted data
% Hsort: the sorted header
% cmp_start: First CMP gather to be corrected
% cmp_end: Last CMP gather to be corrected 
% maxlags: maximum used lag samples
% 
% The output is:
% Dsort: the static corrected seismic data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code if for the book titled: Processing Surface Seismic Reflection
% Data Using MATLAB by Wail A. Mousa & Abdullatif A. Al-Shuhail
% September 2011.

cmp_num=cmp_start:cmp_end;

for i=1:length(cmp_num)	%mincmp and maxcmp are min and max cmp analyzed%
    %Input cmp i
    Da=extracting_cmp(Dsort,Hsort,cmp_num(i));
    [nt,nx]=size(Da);
    Da1=zeros(size(Da));
    Da2=Da1;
    
    strace1=sum(Da,2);	%strace1(i)=first stacked trace of i-th cmp%
    for j=1:nx	%mintrace and maxtrace are min and max traces in i-th cmp%
        ctrace1=my_xcorr(strace1,maxlags,Da(:,j)); %(between -40 and +40 ms)
        [cc,cmax1]=find(ctrace1==max(ctrace1));%pick position of maximum(ctrace1(j))
        Da1(:,j)=[Da(cmax1:nt,j);zeros(cmax1-1,1)];%shift trace(j) by cmax1(j)
      
    end
    
    strace2=sum(Da1,2);	%strace1(i)=first stacked trace of i-th cmp%
    [nt,nx]=size(Da1);
    for j=1:nx	%mintrace and maxtrace are min and max traces in i-th cmp%
        ctrace2=my_xcorr(strace2,maxlags,Da(:,j)); %(between -40 and +40 ms)
        [cc,cmax2]=find(ctrace2==max(ctrace2));%pick position of maximum(ctrace1(j))
        Da2(:,j)=[Da(cmax2:nt,j);zeros(cmax2-1,1)];%shift trace(j) by cmax1(j)
    end
    
    XX=[Hsort.cdp];
    [ii,jj]=find(XX==cmp_num(i));
    Dsort(:,jj)=Da2;
   
end
