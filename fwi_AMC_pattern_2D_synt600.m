% This is modified from fwi_AMC_pattern.m at Jun 17, 2015. The previous one
% is for 3D model building, whereas the previous 2D model for synthetic
% tests is extracted from the 3D AMC model. Now I make it work for 2D
% directly.

path_in = '../../swir_2D_V4_synt600/data/';
% path_in = '../../swir_2D_V4_real680/data/';
file_in = 'swir_30m_2D.vp';

dlen=0.03; NX=2400; NZ=500;
zrange=[0 13.57];
readnew=1; writeout=0;

path_out=path_in;
vpfile=[path_out '/AMC_30m_true_2D.vp'];
vsfile=[path_out '/AMC_30m_true_2D.vs'];

anoamc=-600;

% Plotting option
fontsize=18; labelsize=18;
tsize=21;
fontname='Helvetica';

if readnew == 1    
    fid = fopen ([path_in, file_in], 'r');
    mpin=fread (fid, Inf, 'real*4');
    fclose (fid);
    
    mpin=reshape(mpin,[NZ NX]);
end
xorigin=-36.97-0.3; zorigin=-0.29;
xv=[0:NX-1]*dlen+xorigin; zv=[0:NZ-1]*dlen+zorigin;

zcenamc=220+10;
% xcenamc=1253;
xcenamc=1243; % at 0.3 km for the xorigin=-36.97, also the symmetrical axis in real6[89] series
zupamc=round(0.7/dlen);
zdownamc=round(3/dlen);
xlenamc=round(5/dlen/2);

x1=xcenamc-xlenamc;  x2=xcenamc+xlenamc;

z1=zcenamc-zupamc;  z2=zcenamc+zdownamc;

decayfac=3;
gaufacx=exp(-([-xlenamc:1:xlenamc].^2/(xlenamc.^2)*decayfac));
gaufacy=nan(zupamc+zdownamc+1,1);
gaufacy(1:zupamc+1)=exp(-([-zupamc:1:0].^2/(zupamc.^2)*decayfac));
gaufacy(zupamc+2:end)=exp(-([1:1:zdownamc].^2/(zdownamc.^2)*decayfac));

gaufac=gaufacy*gaufacx;
thres=exp(-decayfac);
gain=1/(1-thres)*anoamc;
for i=1:numel(gaufac)
    if gaufac(i)<=thres
        gaufac(i)=0;
    else
        gaufac(i)=(gaufac(i)-thres)*gain;
    end
end

per=zeros([NZ NX]);
per(z1:z2,x1:x2)=gaufac;

mpout=mpin+per;

if writeout == 1
    fvp=fopen(vpfile,'wb','ieee-le');
    fwrite(fvp,mpout,'real*4',0,'ieee-le');
    fclose(fvp);
    
    [msout,~]=crustal_elastic_parameters_relation(mpout,1);
    fvs=fopen(vsfile,'wb','ieee-le');
    fwrite(fvs,msout,'real*4',0,'ieee-le');
    fclose(fvs);
end

wb=load([path_in '/swir_30m_wb_2D']);
wb=wb*dlen;
wb(:,1)=wb(:,1)+xorigin;
wb(:,2:end)=wb(:,2:end)+zorigin;
wb(:,3)=wb(:,3)+0.3;
shot=load([path_in '/swir_V4_shot_2Dreal_30m']);
shot(:,3)=shot(:,3)*dlen+xorigin;
shot(:,4)=shot(:,4)/1000+zorigin;

figure('Position',[2100 600 900 250],'defaultaxesfontsize',fontsize,...
    'defaultaxesfontname',fontname);

% load ../cmap_hot_jian.mat;
load ../cmap_jet_polar_black_red_central_white.mat
colormap(cmap); caxis([-1 1]*600);

hold on; grid on;box on;
imagesc(xv,zv,squeeze(per(:,:)));
% imagesc(squeeze(per(:,:)));
axis ij tight;
% ylim([0 14]);
daspect([1 1 1]); grid on;
plot(wb(:,1),wb(:,2:3),'k','linewidth',1.);
plot(shot(:,3),shot(:,4)-0.1,'kv','MarkerSize',13,'MarkerFaceColor','g');

colorbar;
set(gca,'layer','top');
set(gca,'xlim',[-1 1]*20,'ylim',zrange);
set(gcf,'paperpositionMode','auto');
xlabel('Across-axis distacne (km) ');
ylabel('Depth (km) ');
set(gca,'XMinorTick','on')
set(gca,'yMinorTick','on')
set(gca,'ytick',[0:4:12],'xtick',[-30:10:100]);
set(gca,'tickdir','out')

load ../cmap_seis_fwi.mat
caxisv=[2500 8300];
% cmap=[1 1 1;cmap];

figure('Position',[1500 500 900 250],'defaultaxesfontsize',fontsize,...
    'defaultaxesfontname',fontname);
imagesc(xv,zv,squeeze(mpin(:,:)));
hold on;
% plot(xcenamc*dlen/1000*[1 1],[0 350],'k','linewidth',2);
colorbar; colormap(cmap);caxis(caxisv);
daspect([1 1 1]); grid on;
plot(wb(:,1),wb(:,2:3),'k','linewidth',1);
plot(shot(:,3),shot(:,4)-0.1,'kv','MarkerSize',13,'MarkerFaceColor','g');
set(gca,'layer','top');
set(gca,'xlim',[-1 1]*20,'ylim',zrange);
set(gcf,'paperpositionMode','auto');
xlabel('Across-axis distance (km) ');
ylabel('Depth (km) ');
set(gca,'XMinorTick','on')
set(gca,'yMinorTick','on')
set(gca,'ytick',[0:4:12],'xtick',[-30:10:100]);
set(gca,'tickdir','out')

figure('Position',[1500 300 900 250],'defaultaxesfontsize',fontsize,...
    'defaultaxesfontname',fontname);
imagesc(xv,zv,squeeze(mpout(:,:)));
hold on;
% plot(xcenamc*dlen/1000*[1 1],[0 350],'k','linewidth',2);
colorbar;  colormap(cmap);caxis(caxisv);
daspect([1 1 1]);  grid on;
plot(wb(:,1),wb(:,2:3),'k','linewidth',1.);
plot(shot(:,3),shot(:,4)-0.1,'kv','MarkerSize',13,'MarkerFaceColor','g');
set(gca,'layer','top');
set(gca,'xlim',[-1 1]*20,'ylim',zrange);
set(gcf,'paperpositionMode','auto');
xlabel('Across-axis distance (km) ');
ylabel('Depth (km) ');
set(gca,'XMinorTick','on')
set(gca,'yMinorTick','on')
set(gca,'ytick',[0:4:12])
set(gca,'xtick',[-30:10:100])
set(gca,'tickdir','out')

figure('Position',[2300 500 300 550]);
plot(squeeze(mpin(:,xcenamc)),zv,'b');
hold on; grid on;
plot(squeeze(mpout(:,xcenamc)),zv,'r');
legend('Raw','AMC','fontname',fontname);
axis ij; grid on
xlim([1000 8000]);
set(gca,'ytick',[0:2:12]);
% colorbar; 
% daspect([1 1 1]); 
set(gca,'fontsize',fontsize,'fontname',fontname);
set(gca,'layer','top');
set(gcf,'paperpositionMode','auto');
