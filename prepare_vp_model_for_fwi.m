% we need segy files+ pick files for every segy file+ velocity file with
% xyz format.
%we need to change the range of velocity model on figure axis
%picked in line 228 data sould be chanaged based on the max shot
clear all
close all
clc
cd ../
cd ('/Volumes/DRIVE_DATA/TEST3/')
obs_number = 2
shot_depth_meter = 14
max_num_shot = 182
%________________________________________
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
%%to run this part we need out.dws file, but it is not important!
% figure1 = figure
% set(0, 'defaultfigurecolor', [1 1 1])
%  set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 25, 13], ...
%      'PaperUnits', 'centimeters', 'PaperSize', [18, 9])
% pcolor(Distance, fliplr(Depth), fliplr(Velocity_raw));
% colormap(flipud(jet));
% fs=20;
% set(gca,'FontName','Times New Roman','fontsize', fs)
% shading interp
% set(gca, 'YDir', 'reverse')
% axis equal
% image=colormap(fliplr(jet));
% set(gca,'xaxislocation', 'top')
% set(gca,'YDir','reverse')
% set(gca,'FontName','Times New Roman','fontsize',fs)
% xlabel('Distance [km]','FontName','Times New Roman','FontSize',fs); ylabel('Depth [km]','FontName','Times New Roman','FontSize',fs);
% axis([0 43 0, 15])
% caxis([1.6 8.0]);
% set(gca,'TickDir','out'); 
% image(1,:) = [1 1 1];
% colormap(image);
% 
% cut=importdata('out.dws');
% I=min(cut(:,1));
% J=1;
% dws=[];
% while I<=max(cut(:,1))
%     II=find(cut(:,1)==I);
%     cut_distance(J)=I;
%     cut_depth(1:length(II),J)=sort(cut(II(1):II(end),2));
%     cut_ave(1:length(II),J)=(cut(II(1):II(end),3));
%     JJ=find(cut(:,1)>I,1);
%     I=cut(JJ,1);
%     J=J+1;
% end
% 
% for I=1:length(cut_ave(1,:))
%     L2=find(cut_ave(:,I)~=0);
%     if L2>0
%         Q(I)=cut_depth(L2(end),I);
%     else
%         Q(I)=cut_depth(1,I);
%     end
% end
% 
% hold on
% plot(cut_distance+2.1,fliplr(Q),'w--', 'linewidth', 2)
% Off=meshgrid(Distance);
% Off=Off(1:length(Depth(:,1)),:);
% [C, cfig]=contourm(fliplr(Depth),Off,fliplr(Velocity_raw), [2.5  3.5  4.5  5.5  6.5  7.5 8], 'Color', [1 1 1],'LineWidth', 1.5, 'ShowText', 'on');
% ht=clabelm(C, cfig);
% set(ht,'FontSize',fs-10);
% set(gca,'FontName','Times New Roman','fontsize',fs)
% c=colorbar;set(gca,'FontName','Times New Roman','fontsize',fs)
% ylabel(c,'V_P [km/s]','FontName','Times New Roman','FontSize', fs)
% hold on
% set(gcf, 'PaperPositionMode', 'auto' );
% lab=strcat('2D OBS New Velocity with Contours_AVE');
% print(gcf,lab,'-dpng','-r300')
% saveas(figure1,'fig1.png')
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


close all
figure2 = figure
set(0, 'defaultfigurecolor', [1 1 1])
 set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 25, 13], ...
     'PaperUnits', 'centimeters', 'PaperSize', [18, 9])
% pcolor(Distance, fliplr(Depth), fliplr(Velocity_raw));
pcolor(Distance, (Depth), (Velocity_raw));
colormap(flipud(jet));
hCB=colorbar;
hCB.Title.String='Vp (km/s)';
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
axis([0 160 0, 15])
caxis([1.6 8.0]);
set(gca,'TickDir','out'); 
image(1,:) = [1 1 1];
colormap(image);
% c=colorbar;set(gca,'FontName','Times New Roman','fontsize',fs)
% ylabel(c,'V_P [km/s]','FontName','Times New Roman','FontSize', fs)
title('Original')
%D gives the number of points in grid (horizontally and vertically)
D=xq.'; 
vel_new=vel_new.';
%% Extract water bottom file
wat_bot_grid=zeros(1,length(vel_new(1,:)));
wat_bot_depth=wat_bot_grid;
for I=1:length(vel_new(1,:))
    K=find(vel_new(:,I)>1.6,1);
    wat_bot_grid(I)=K;
    wat_bot_depth(I)=D(K,I);
end
saveas(figure2,'fig2.png')

%% QC
figure3 = figure
set(0, 'defaultfigurecolor', [1 1 1])
 set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 25, 13], ...
     'PaperUnits', 'centimeters', 'PaperSize', [18, 9])
% pcolor(dist, fliplr(D), fliplr(vel_new));
pcolor(dist, (D), (vel_new));
colormap(flipud(jet));
hCB=colorbar;
hCB.Title.String='Vp (km/s)';
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
axis([0 160 0, 15])
caxis([1.6 8.0]);
set(gca,'TickDir','out'); 
image(1,:) = [1 1 1];
colormap(image);
% c=colorbar;set(gca,'FontName','Times New Roman','fontsize',fs)
% ylabel(c,'V_P [km/s]','FontName','Times New Roman','FontSize', fs)
title('Resampled')
hold on
plot(dist, (wat_bot_depth),'k', 'linewidth',2)
saveas(figure3,'fig3.png')
wat_depth=[(1:length(wat_bot_grid)).', (wat_bot_grid).'];

savefile = 'eswir_water_bottom_grid';
fid=fopen(savefile,'wt');
for ll=1:length(wat_bot_grid)
    fprintf(fid, '%5d\t%4d\n',wat_depth(ll,:));
end
    fclose(fid);    
% dlmwrite('eswir_water_bottom_grid', wat_depth, 'precision', 6, 'delimiter', '\t')
% type eswir_water_bottom_grid
system('mv /Volumes/DRIVE_DATA/TEST3/eswir_water_bottom_grid /Volumes/DRIVE_DATA/TEST3/Preprocessed_files')
%% Receiver location file in grid points and depth in meters

format long g

A=importdata('picks.2.dat');
shot_loc=A(2:end,1);
gp=zeros(size(shot_loc));
for I=1:length(shot_loc)
    J=find(dist>=shot_loc(I),1);
    gp(I)=J;
end
rec_depth=[(1:max_num_shot).', gp,(shot_depth_meter*ones(size(1:max_num_shot))).'];

lab1=strcat('eswir_rec_depth_file');
savefile = lab1;
fid=fopen(savefile,'wt');
for ll=1:max_num_shot
    fprintf(fid, '%7d\t%7d\t%8.1f\n',rec_depth(ll,:));
end
fclose(fid);
system('mv /Volumes/DRIVE_DATA/TEST3/eswir_rec_depth_file /Volumes/DRIVE_DATA/TEST3/Preprocessed_files')
 
% dlmwrite('eswir_rec_depth_file', rec_depth, 'precision', 6, 'delimiter', '\t')
% type eswir_rec_depth_file
% clear shot_num
%% Source depth file
shot_num=[];
for J=1:obs_number
    lab=strcat('picks.', num2str(J),'.dat');
    A=importdata(lab);
    K=find(dist>=A(1,1),1);
    I=find(D(:,K)>=A(1,2),1);
    shot_num(J,1)=J;
    shot_num(J,2)=J;
    shot_num(J,3)=K;
    shot_num(J,4)=A(1,2)*1000;
    shot_num(J,5)=1;
end
sou_depth=[shot_num];
dlmwrite('/Volumes/DRIVE_DATA/TEST3/Preprocessed_files/eswir_shot_depth_file', shot_num, 'precision', 6, 'delimiter', '\t')
%type eswir_shot_depth_file

s=scatter((dist(shot_num(:,3))),((shot_num(:,4)))./1000, 90,'MarkerEdgeColor',[0 0 0],...
    'MarkerFaceColor',[1 1 1],'LineWidth',2);

s=scatter(dist(rec_depth(1:10:end,2)),rec_depth(1:10:end,3)./1000,'*','r',...
'LineWidth',0.2);
%% Prepare Vp files

SF=zeros(size(vel_new));
rho=SF;
vel_shear=SF;
for I=1:length(vel_new(1,:))
    vel_shear(:,I)=vel_new(:,I)./sqrt(3);
    II=find(vel_new(:,I)>1.6,1);
    vel_shear(1:II-1,I)=0;
    for J=1:length(SF(:,I))
        if vel_new(J,I)>7.95
            SF(J,I)=0;
            rho(J,I)=3.3;
        elseif vel_new(J,I)<5.07
            SF(J,I)=1;
            rho(J,I)=1;
        elseif vel_new(J,I)==5.07
            rho(J,I)=2.51;
        else
            SF(J,I)=((vel_new(J,I)-7.95)./-2.88);
            rho(J,I)=3.3-(0.79*SF(J,I));
        end
    end
end

vel_new=vel_new.*1000;
vel_shear=vel_shear.*1000;

vp=zeros(length(vel_new(:,1))*length(vel_new(1,:)),3);
vs=zeros(length(vel_new(:,1))*length(vel_new(1,:)),3);
dens=zeros(length(vel_new(:,1))*length(vel_new(1,:)),3);

Offs=zeros(size(D));
for I=1:length(Offs(:,1))
    Offs(I,:)=dist;
end
vp=[Offs(:).*1000, D(:).*1000, vel_new(:)];
vs=[Offs(:).*1000, D(:).*1000, vel_shear(:)];
dens=[Offs(:).*1000, D(:).*1000, rho(:)];

dlmwrite('/Volumes/DRIVE_DATA/TEST3/Preprocessed_files/eswir_vp_input.vp', vp, 'precision', 6, 'delimiter', '\t')
% type eswir_vp_input

id=fopen(lab, 'w');
fwrite(id,A,'real*4');
fclose(id);

fileID = fopen('/Volumes/DRIVE_DATA/TEST3/Preprocessed_files/eswir_vp_input','wb','ieee-le');
fwrite(fileID, vp(:,3),'real*8',0,'ieee-le');
fclose(fileID);

dlmwrite('/Volumes/DRIVE_DATA/TEST3/Preprocessed_files/eswir_vs_input.vs', vs, 'precision', 6, 'delimiter', '\t')
% type eswir_vs_input

fileID = fopen('/Volumes/DRIVE_DATA/TEST3/Preprocessed_files/eswir_vs_input','w');
fwrite(fileID, vs(:,3),'real*4');
fclose(fileID);

% dlmwrite('eswir_rho_input.rho', rho, 'precision', 6, 'delimiter', '\t')
% fileID = fopen('eswir_rho_input','w');
% fwrite(fileID, dens(:,3),'single');
% fclose(fileID);

fileID = fopen('/Volumes/DRIVE_DATA/TEST3/Preprocessed_files/eswir_rho_input','w');
fwrite(fileID, dens(:,3).*1000,'real*4');
fclose(fileID);
% type eswir_rho_input
%% Define mute file
for nn=1:obs_number
    lab=strcat('0',num2str(nn),'b_obs',num2str(nn),'.segy');
    SEISMIC=read_segy_file(lab);
    data=SEISMIC.traces(:,1:length(SEISMIC.traces(1,:)));
    data_cut=data;
    mut_top=[];
    mut_bot=[];
    for I=1:length(data(1,:))
        J=find(data(:,I)~=0,1);
        if isempty(J)~=1
            mut_top(I)= J-25;
        else
            mut_top(I)=mut_top(I-1);
        end
        ind=round(abs(SEISMIC.headers(4,(I)))/35);
        mut_bot(I)=mut_top(I)+ind;
        data_cut(mut_bot(I):end,I)=0;
        data_cut(1:mut_top(I),I)=0;
    end
    nnn=find(SEISMIC.headers(4,:)>=0,1);
    mut_bot(nnn-15:nnn+20)=mut_top(nnn-15:nnn+20)-95;
    close all
    imagesc(data);
    caxis([-1e5 1e5])
    colormap gray
    hold on
    plot(1:max_num_shot,mut_top,'r')
    plot(1:max_num_shot,mut_bot,'r')
%     figure was ok
    mut_samp=[(1:max_num_shot).', mut_top.' ];
    lab1=strcat('/Volumes/DRIVE_DATA/TEST3/Preprocessed_files/eswi_tmut_', num2str(nn),'hydr');
    dlmwrite(lab1, mut_samp, 'precision', 6, 'delimiter', '\t')
    type(lab1)
    clear mut_samp
    
    mut_samp=[(1:max_num_shot).', mut_bot.' ];
    lab1=strcat('/Volumes/DRIVE_DATA/TEST3/Preprocessed_files/eswi_bmut_', num2str(nn),'hydr');
    dlmwrite(lab1, mut_samp, 'precision', 6, 'delimiter', '\t')
    type(lab1)
    clear mut_samp
    mut_samp=[(1:max_num_shot).', 2*(mut_top.'), 2*(mut_bot).'];
    lab1=strcat('/Volumes/DRIVE_DATA/TEST3/Preprocessed_files/eswi_mut_000', num2str(nn),'hydr');
    if nn>=10
        lab1=strcat('/Volumes/DRIVE_DATA/TEST3/Preprocessed_files/eswi_mut_00', num2str(nn),'hydr');
    end
    savefile = lab1;
    fid=fopen(savefile,'wt');   
    for ll=1:max_num_shot
         fprintf(fid, '%7d\t%7d\t%7d\n',mut_samp(ll,:));
    end    
    fclose(fid);    
%   fprintf(fid,'%3d\n',mut_samp)
%   dlmwrite(lab1, mut_samp, 'precision', 6, 'delimiter', '\t')
    type(lab1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SEISMIC.traces=data_cut;
    SEISMIC=s_resample(SEISMIC,2);
    
    off=SEISMIC.headers(4,:)./1000;
    t=0:(SEISMIC.step)./1000:(SEISMIC.last)./1000;
    
    close all
    fs=16;
    for mm=1:length(SEISMIC.traces(1,:))
        SEISMIC.traces(1:mut_top(mm)*2,mm)=0;
        SEISMIC.traces(mut_bot(mm)*2:end,mm)=0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SEISMIC=s_filter(SEISMIC,{'ormsby',[1 3 10 25]});
    %figure
    hold on
    set(0,'defaultfigurecolor',[1 1 1])
    t=(0:SEISMIC.step/1000:SEISMIC.last/1000);
    mwigb(fliplr(SEISMIC.traces),1.0,-1*fliplr(off),t)
%     mwigb(fliplr(data_cut),1.0,-1*fliplr(off),t)
    set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 24, 15], ...
        'PaperUnits', 'centimeters', 'PaperSize', [12, 9])
    set(gca,'TickDir','out');
    set(gca,'xaxislocation','top');
    set(gca,'TickDir','out');
    fs=16;
    ylim([2.5 12])
    hold on
    plot(-1*(off),SEISMIC.step/1000* 2*(mut_top),'r')
    plot(-1*(off),SEISMIC.step/1000* 2*(mut_bot),'r')
    %figure was not ok
    set(gca,'fontsize', fs)
    set(gca,'FontName','Times New Roman')
    xlabel('Signed Offset [km]', 'FontName','Times New Roman','FontSize', fs)
    ylabel('Time [s]', 'FontName','Times New Roman','FontSize', fs)
    lab=strcat('Fig obs', num2str(nn),'_input_with_mutes');
    set(gcf, 'PaperPositionMode', 'auto' );
    print(gcf,'-r300','-dpng',lab)
    
    KK=find(all(ismember(SEISMIC.header_info,'ffid'),4)==1);
    SEISMIC.headers(KK,:)=nn;
    lab1=strcat('eswir_ffid_',num2str(nn),'hydr.segy');
    write_segy_file(SEISMIC,lab1);
end
%system('mv /Volumes/DRIVE_DATA/TEST3/eswi_mut_* /Volumes/DRIVE_DATA/TEST3/Preprocessed_files')
%ADD by yaser
system('mv /Volumes/DRIVE_DATA/TEST3/eswir_ffid_* /Volumes/DRIVE_DATA/TEST3/Preprocessed_files')

%% Convert segy to binary
for nn=1:obs_number
    %yaser changed the path
    lab=strcat('/Volumes/DRIVE_DATA/TEST3/Preprocessed_files/eswir_ffid_',num2str(nn),'hydr.segy');
    SEISMIC=read_segy_file(lab);
    data=SEISMIC.traces(:,1:length(SEISMIC.traces(1,:)));   
    if nn<10
        lab2=strcat('eswi_ffid_000',num2str(nn),'hydr');
    else
        lab2=strcat('eswi_ffid_00',num2str(nn),'hydr');
    end 
    fileID = fopen(lab2,'w');
    fwrite(fileID, data,'real*4');
    fclose(fileID);
end
system('mv /Volumes/DRIVE_DATA/TEST3/eswi_ffid_* /Volumes/DRIVE_DATA/TEST3/Preprocessed_files')

fid=fopen('/Volumes/DRIVE_DATA/TEST3/Preprocessed_files/eswi_ffid_0001hydr','r')
A=fread(fid, Inf, 'float');
B=reshape(A, [length(data(:,1)), length(data(1,:))]);
fclose(fid);
%% Put shot gathers in 1 file

lab=strcat('/Volumes/DRIVE_DATA/TEST3/Preprocessed_files/eswir_ffid_1hydr.segy');
SEISMIC=read_segy_file(lab);
A=SEISMIC;

data=zeros(length(SEISMIC.traces(:,1)),10*length(SEISMIC.traces(1,:)));
headers=zeros(length(SEISMIC.headers(:,1)),10*length(SEISMIC.traces(1,:)));
for I=1:obs_number   
%     if nn<10
        lab2=strcat('/Volumes/DRIVE_DATA/TEST3/Preprocessed_files/eswir_ffid_',num2str(I),'hydr');
%     else
%         lab2=strcat('eswi_ffid_00',num2str(nn),'hydr');
%     end
    SEISMIC=read_segy_file(lab);
    J=length(SEISMIC.traces(1,:));
    data(:,(I*J)-J+1:I*J)=SEISMIC.traces;
    headers(:,(I*J)-J+1:I*J)=SEISMIC.headers;
    headers(2,(I*J)-J+1:I*J)=I;
end
A.traces=zeros(size(data));
A.traces=data;
A.headers=zeros(size(headers));
A.headers=headers;
filename=('/Volumes/DRIVE_DATA/TEST3/Preprocessed_files/eswir_ffid_1_10.segy');
write_segy_file(A,filename); 