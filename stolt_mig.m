function Dmigrated=stolt_mig(Dstacked,v,dt,dx,t)
% This code is written to perform Stolt time migration on a stacked seismic data. 
% The required inputs is:
% Dstacked: the stacked seismic data
% v: the velocity to be used in the migration process
% dt: time sampling interval
% dx: spatial sampling interval
% t: a vector contianing the time values 
% The output is:
% Dmigrated: the time migrated data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code if for the book titled: Processing Surface Seismic Reflection
% Data Using MATLAB by Wail A. Mousa & Abdullatif A. Al-Shuhail
% September 2011.


[nt,ncdp]=size(Dstacked);
num_f_pts=nt;
num_pts=num_f_pts;

tz=t*v;
tz=repmat(tz,ncdp,1)';
Dstacked=Dstacked.*tz;

kx=(2*pi*linspace(-0.5,0.5,num_pts)/dx)';
U_w_kx=fft2(Dstacked,num_f_pts,num_pts);
U_kz_kx=zeros(size(U_w_kx));

for i=1:num_pts
    U_kz_kx(:,i)=U_w_kx(:,i).*sqrt(1-((v*kx).^2)./(2*U_w_kx(:,i)));
    U_kz_kx(:,i)=sin(U_kz_kx(:,i));
    S=(v/2)*U_kz_kx(:,i)./sqrt(kx.^2+U_kz_kx(:,i).^2);
    U_kz_kx(:,i)=U_kz_kx(:,i).*S; 
end
Dmigrated=real(ifft2(U_kz_kx));
Dmigrated=Dmigrated(:,1:ncdp);
Dmigrated=AGCgain(Dmigrated,dt,0.5,2);