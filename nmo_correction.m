function [Dsort,Hsort]=nmo_correction(Dsort,Hsort,vnmo,tnmo,cmp_start,cmp_end,cmp_step)
% This code is written to perform NMO-correction for CMP seismic gathers.
% The required inputs are:
% Dsort: the sorted data
% Hsort: the sorted header
% vnmo: NMO picked velocities
% tnmo: NMO picked time values
% cmp_start: First CMP gather to be corrected
% cmp_end: Last CMP gather to be corrected 
% cmp_step: CMP gather step size to be used
% 
% The output is:
% Dsort: the NMO-corrected seismic data
% Hsort: the NMO-corrected header (nothing changed just to re-store it)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code if for the book titled: Processing Surface Seismic Reflection
% Data Using MATLAB by Wail A. Mousa & Abdullatif A. Al-Shuhail
% September 2011.

cmp_num=cmp_start:cmp_step:cmp_end;
for kkk=1:length(cmp_num)
    Max_Stretch = 50;
    [Da,dt,t,cdp,jj,cmp_offset1]=extracting_cmp(Dsort,Hsort,cmp_num(kkk));
    h =cmp_offset1;
    % NMO
    Da = nmo(Da,dt,h,tnmo(kkk,:),vnmo(kkk,:),Max_Stretch);
    Dsort(:,jj)=Da;
    if cmp_num(kkk)<400
        for ll=1:cmp_step-1
            [Da,dt,t,cdp,jj,cmp_offset1]=extracting_cmp(Dsort,Hsort,cmp_num(kkk)+ll);
            disp(['I am CMP number: ',num2str(cmp_num(kkk)+ll),''])
            h =cmp_offset1;
            Da = nmo(Da,dt,h,tnmo(kkk,:),vnmo(kkk,:),Max_Stretch);
            Dsort(:,jj)=Da;
        end
    end
end