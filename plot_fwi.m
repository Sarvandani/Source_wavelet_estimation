clear all
close all
clc
format long g
cd ../
cd D:\KEEZ\PhD_Momoh_MC\Tomography\OBS_Data\Data_Preparation_for_FWI\Ek

fid=fopen('winv.mp00', 'r');
v=fread(fid, Inf,'float');
fclose(fid);
vv=reshape(v,[667,1467]);
figure
imagesc(vv)
caxis([1500 8000])
title('mp00')

fid=fopen('winv.mp01', 'r');
v=fread(fid, Inf,'float');
fclose(fid);
vv=reshape(v,[667,1467]);
figure
imagesc(vv)
caxis([1500 8000])
title('mp01')

fid=fopen('winv.mp02', 'r');
v=fread(fid, Inf,'float');
fclose(fid);
vv=reshape(v,[667,1467]);
figure
imagesc(vv)
caxis([1500 8000])
title('mp02')

fid=fopen('winv.mp03', 'r');
v=fread(fid, Inf,'float');
fclose(fid);
vv=reshape(v,[667,1467]);
figure
imagesc(vv)
caxis([1500 8000])
title('mp03')

% fid=fopen('winv.mp04', 'r');
% v=fread(fid, Inf,'float');
% fclose(fid);
% vv=reshape(v,[667,1467]);
% figure
% imagesc(vv)
% caxis([1500 8000])
% title('mp04')


%%
% A=read_segy_file('eswir_ffid_5hydr.segy');
% Distance=A.headers(4,:)/1000;
% time=linspace(0,16,8001);
sca=zeros([10,2]);
for I=1:10
    if I<10
        lab=strcat('dmod.hydr.s00', num2str(I),'.ite001');
    elseif I>=10
        lab=strcat('dmod.hydr.s0', num2str(I),'.ite001');
    end
    lab2=strcat('obs', num2str(I), '_decon_filtered_resampled.segy');
    A=read_segy_file(lab2);
    Distance=A.headers(4,:)/1000;
    time=linspace(0,16,8001);
    fid=fopen(lab, 'r');
    v=fread(fid,inf,'float');
    fclose(fid);
    vv=reshape(v,[182,8001]).';
    scaler=max(abs(A.traces(:))/max(abs(vv(:))));
    sca(I,2)=scaler;
    sca(I,1)=I;
    figure
    set(0,'defaultfigurecolor',[1 1 1])
    mwigb((vv),1,(Distance(1,:))/1,time)
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
    title('output')
    lab=strcat('data_mod', num2str(I));
    print(gcf,'-r300','-dpng',lab)
    figure
    set(0,'defaultfigurecolor',[1 1 1])
    mwigb((A.traces),1,(Distance(1,:))/1,time)
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
    title('input')
    lab=strcat('data_real', num2str(I));
    print(gcf,'-r300','-dpng',lab)
end
savefile = 'shot_norm_file';
fid=fopen(savefile,'wt');
for ll=1:10
    fprintf(fid, '%5d\t%10f\n',sca(ll,:));
end
fclose(fid);
type shot_norm_file
system('copy D:\KEEZ\PhD_Momoh_MC\Tomography\OBS_Data\Data_Preparation_for_FWI\Ek\shot_norm_file D:\KEEZ\PhD_Momoh_MC\Tomography\OBS_Data\Data_Preparation_for_FWI\Preprocessed_files')

% fid=fopen('eswi_ffid_0005hydr', 'r');
% v=fread(fid,inf,'float');
% fclose(fid);
% vv2=reshape(v,[182,8001]);
figure
set(0,'defaultfigurecolor',[1 1 1])
mwigb(fliplr(A.traces),1,(Distance(1,:))/1,time)
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
title('input')


figure
hold on
plot(A.traces(:,110))
plot(vv(:,110), 'r')
plot(vv(:,110)*-scaler, 'g')
% dlmwrite('eswir_water_bottom_grid', wat_depth, 'precision', 6, 'delimiter', '\t')
% type eswir_water_bottom_grid
% for I=5
%     if I<10
%         lab=strcat('dres.mplx.hydr.s00', num2str(I));
% %     elseif I>=10
% %         lab=strcat('dmod.mplx.hydr.s0', num2str(I));
%     end
%     fid=fopen(lab, 'r');
%     v=fread(fid,inf,'float');
%     fclose(fid);
%     vv3=reshape(v,[182,8001]).';
%     figure
%     set(0,'defaultfigurecolor',[1 1 1])
%     mwigb(fliplr(vv3),1,(Distance(1,:))/1,time)
%     set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 24, 15], ...
%         'PaperUnits', 'centimeters', 'PaperSize', [18, 9])
%     set(gca,'TickDir','out');
%     set(gca,'xaxislocation','top');
%     set(gca,'TickDir','out');
%     fs=16;
%     ylim([2.5 12])
%     set(gca,'fontsize', fs)
%     set(gca,'FontName','Times New Roman')
%     xlabel('Trace number', 'FontName','Times New Roman','FontSize', fs)
%     ylabel('Time [s]', 'FontName','Times New Roman','FontSize', fs)
%     set(gcf, 'PaperPositionMode', 'auto' );
%     labs=strcat('dmod');
%     title('output')
% %     print(gcf,'-r300','-dpng',labs)
% end