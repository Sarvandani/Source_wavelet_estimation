clear all
close all
clc
cd ../
cd D:\KEEZ\PhD_Momoh_MC\Tomography\OBS_Data\Data_Preparation_for_FWI
vel_ave=importdata('ave.xyz');
% vel_raw=importdata('smesh_sorted.xyz');
vel_raw=importdata('ave.xyz');
Offset=vel_ave(:,1);
Depth=[];
Distance=[];
Velocity_ave=[];
Velocity_raw=[];
I=min(vel_ave(:,1));
J=1;
while I<=max(vel_ave(:,1))
    II=find(vel_ave(:,1)==I);
    Distance(J)=I;
    Depth(1:length(II),J)=sort(vel_ave(II(1):II(end),2));
    K=find(vel_ave(II,2)==Depth(1:length(II),J));
    Velocity_ave(1:length(II),J)=sort(vel_ave(II(1):II(end),3));
    K=find(Velocity_ave(1:length(II),J)==max(Velocity_ave(1:length(II),J)));
    Velocity_ave(K(1):end,J)=max(Velocity_ave(1:length(II),J));
    mut(J)=(K(1));
    JJ=find(vel_ave(:,1)>I,1);
    I=vel_ave(JJ,1);
    J=J+1;
end

for I=1:length(Depth(1,:))
     II=find(Depth(:,I)==max(Depth(:,I)));
    Depth(II(1):end,I)=max(Depth(:,I))+(0:0.125:0.125*(length(Depth(II(1):end,I))-1));
end
I=min(vel_raw(:,1));
J=1;
while I<=max(vel_raw(:,1))
    II=find(vel_raw(:,1)==I);
    Distance(J)=I;
    Depth(1:length(II),J)=sort(vel_raw(II(1):II(end),2));
    Velocity_raw(1:length(II),J)=sort(vel_raw(II(1):II(end),3));
    K=find(Velocity_raw(1:length(II),J)==max(Velocity_raw(1:length(II),J)));
    Velocity_raw(K(1):end,J)=max(Velocity_raw(1:length(II),J));
    mut(J)=(K(1));
    JJ=find(vel_raw(:,1)>I,1);
    I=vel_raw(JJ,1);
    J=J+1;
end
for I=1:length(Depth(1,:))
     II=find(Depth(:,I)==max(Depth(:,I)));
     Depth(II(1):end,I)=max(Depth(:,I))+(0:0.125:0.125*(length(Depth(II(1):end,I))-1));
end

for I=1:length(Velocity_raw(1,:))
J=find(Velocity_raw(:,I)==max(Velocity_raw(:,I)),1);
Velocity_raw(J:end,I)=max(Velocity_raw(:,I));
end

figure
set(0, 'defaultfigurecolor', [1 1 1])
 set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 25, 10], ...
     'PaperUnits', 'centimeters', 'PaperSize', [18, 9])
pcolor(Distance, fliplr(Depth), fliplr(Velocity_raw));
colormap(flipud(jet));
fs=20;
set(gca,'FontName','Times New Roman','fontsize', fs)
shading interp
set(gca, 'YDir', 'reverse')
axis equal
image=colormap(fliplr(jet));
set(gca,'xaxislocation', 'top')
set(gca,'YDir','reverse')
set(gca,'FontName','Times New Roman','fontsize',fs)
xlabel('Distance [km]','FontName','Times New Roman','FontSize',fs); ylabel('Depth [km]','FontName','Times New Roman','FontSize',fs);
axis([0 43 0, 15])
caxis([1.6 8.0]);
set(gca,'TickDir','out'); 
image(1,:) = [1 1 1];
colormap(image);

cut=importdata('out.dws');
I=min(cut(:,1));
J=1;
dws=[];
while I<=max(cut(:,1))
    II=find(cut(:,1)==I);
    cut_distance(J)=I;
    cut_depth(1:length(II),J)=sort(cut(II(1):II(end),2));
    cut_ave(1:length(II),J)=(cut(II(1):II(end),3));
    JJ=find(cut(:,1)>I,1);
    I=cut(JJ,1);
    J=J+1;
end

for I=1:length(cut_ave(1,:))
L2=find(cut_ave(:,I)~=0);
if L2>0
Q(I)=cut_depth(L2(end),I);
else
    Q(I)=cut_depth(1,I);
end
end

hold on
plot(cut_distance+2.1,fliplr(Q),'w--', 'linewidth', 2)
Off=meshgrid(Distance);
Off=Off(1:length(Depth(:,1)),:);
[C, cfig]=contourm(fliplr(Depth),Off,fliplr(Velocity_raw), [2.5  3.5  4.5  5.5  6.5  7.5 8], 'Color', [1 1 1],'LineWidth', 1.5, 'ShowText', 'on');
ht=clabelm(C, cfig);
set(ht,'FontSize',fs-10);
set(gca,'FontName','Times New Roman','fontsize',fs)
c=colorbar;set(gca,'FontName','Times New Roman','fontsize',fs)
ylabel(c,'V_P [km/s]','FontName','Times New Roman','FontSize', fs)
hold on
% plot(flipud(Offset1)+0.995,(gravity_crust)*-1./1000, 'r--', 'linewidth',4);
set(gcf, 'PaperPositionMode', 'auto' );
lab=strcat('2D OBS New Velocity with Contours_AVE');
% print(gcf,lab,'-dpng','-r300')
%% Interpolate data

[I,J]=find(Depth<0);
new_depth=Depth((I(end)+1):end,:);
new_velocity=Velocity_raw((I(end)+1):end,:);

dist=Distance(1):30/1000:Distance(end);
dep=(new_depth(1,1):30/1000:new_depth(end,1));
d=meshgrid(Distance);
d=d(1:length(new_depth(:,1)),:);

[xq,yq]=meshgrid(dep,dist);
vel_new = griddata(new_depth,d,new_velocity,xq,yq);
D=xq.'; 
vel_new=vel_new.';
%% QC
close all
figure
set(0, 'defaultfigurecolor', [1 1 1])
 set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 25, 10], ...
     'PaperUnits', 'centimeters', 'PaperSize', [18, 9])
% pcolor(dist, fliplr(D), fliplr(vel_new));
pcolor(dist, (D), (vel_new));
colormap(flipud(jet));
fs=20;
set(gca,'FontName','Times New Roman','fontsize', fs)
shading interp
set(gca, 'YDir', 'reverse')
axis equal
image=colormap(fliplr(jet));
set(gca,'xaxislocation', 'top')
set(gca,'YDir','reverse')
set(gca,'FontName','Times New Roman','fontsize',fs)
xlabel('Distance [km]','FontName','Times New Roman','FontSize',fs); ylabel('Depth [km]','FontName','Times New Roman','FontSize',fs);
axis([0 43 0, 15])
caxis([1.6 8.0]);
set(gca,'TickDir','out'); 
image(1,:) = [1 1 1];
colormap(image);
title('Resampled')

fid=fopen('winv.mp01', 'r');
v=fread(fid, Inf,'float');
fclose(fid);
vv=reshape(v,[667,1467]);
figure
set(0, 'defaultfigurecolor', [1 1 1])
 set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 25, 10], ...
     'PaperUnits', 'centimeters', 'PaperSize', [18, 9])
% pcolor(dist, fliplr(D), fliplr(vel_new));
pcolor(dist, (D), (vv)./1000);
colormap(flipud(jet));
fs=20;
set(gca,'FontName','Times New Roman','fontsize', fs)
shading interp
set(gca, 'YDir', 'reverse')
axis equal
image=colormap(fliplr(jet));
set(gca,'xaxislocation', 'top')
set(gca,'YDir','reverse')
set(gca,'FontName','Times New Roman','fontsize',fs)
xlabel('Distance [km]','FontName','Times New Roman','FontSize',fs); ylabel('Depth [km]','FontName','Times New Roman','FontSize',fs);
axis([0 43 0, 15])
caxis([1.6 8.0]);
set(gca,'TickDir','out'); 
image(1,:) = [1 1 1];
colormap(image);
title('Waveform inversion')

time=linspace(0,16,8001);
Distance=1:1:182;
% fid=fopen('dmod.hydr.s001.ite001', 'r');
fid=fopen('dmod.mplx.hydr.s001', 'r');
v=fread(fid,inf,'float');
fclose(fid);
vv=reshape(v,[182,8001]).';
figure
set(0,'defaultfigurecolor',[1 1 1])
mwigb(fliplr(vv),1.0,(Distance(1,:))/1,time)
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 24, 15], ...
    'PaperUnits', 'centimeters', 'PaperSize', [18, 9])
set(gca,'TickDir','out');
set(gca,'xaxislocation','top');
set(gca,'TickDir','out');
fs=16;
ylim([2.5 12])
set(gca,'fontsize', fs)
set(gca,'FontName','Times New Roman')
xlabel('Trace number', 'FontName','Times New Roman','FontSize', fs)
ylabel('Time [s]', 'FontName','Times New Roman','FontSize', fs)
set(gcf, 'PaperPositionMode', 'auto' );
% print(gcf,'-r300','-dpng',lab)

v=read_segy_file('eswir_ffid_5hydr.segy');
v=v.traces;
figure
set(0,'defaultfigurecolor',[1 1 1])
mwigb(fliplr(vv),1.0,(Distance(1,:))/1,time)
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 24, 15], ...
    'PaperUnits', 'centimeters', 'PaperSize', [18, 9])
set(gca,'TickDir','out');
set(gca,'xaxislocation','top');
set(gca,'TickDir','out');
fs=16;
ylim([2.5 12])
set(gca,'fontsize', fs)
set(gca,'FontName','Times New Roman')
xlabel('Trace number', 'FontName','Times New Roman','FontSize', fs)
ylabel('Time [s]', 'FontName','Times New Roman','FontSize', fs)
set(gcf, 'PaperPositionMode', 'auto' );
% print(gcf,'-r300','-dpng',lab)