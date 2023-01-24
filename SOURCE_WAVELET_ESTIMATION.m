%% I would like to acknowledge I would like to acknowledge the authors and team of several used functions especially SeisLab, 
%% ,Wail A. Mousa & Abdullatif A. Al-Shuhail 2011, E. Rietsch and M.D.Sacchi.   
%% Author: Sarvandani
clear all
clc
%% this parameter is the maximum range of sampled time for decimating or exported file:
max_sample_time = 1001;
%% scaling factor for source
SCALE_FACTOR_first_instrument = 10.^2;
% this X1 is offset lines or column in the data set. You
% can check it by header(:,:) = seismic.headers.
X1 = 1;
%% this change the range of amplitude and make different contrast on the figure
contrast_coff = 10000;
%% this is the maximum range of time on signals in terms of sec for autocorrelations
max_lag = 2;
Number_of_segy_files = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reading data
for a=1:Number_of_segy_files
     lab=strcat('synthetic_data',num2str(a),'.segy');
     seismic(a)=read_segy_file(lab);
%      seismic=s_resample(seismic,resampling_cofficient);
     dt = seismic.step;
   %% dt must be in terms of second
     dt = dt/1000;
     off1=seismic(a).headers(X1,:)./1000;
     t=(0:seismic(a).step/1000:seismic(a).last/1000);    
   
end
%%%%%%%%%%%
%% all_offset1 repeats the first iteration two times and we reomve the extra 
%% in all_offset_all
all_offset1 = off1;
for aa = 1:Number_of_segy_files
    offsetToAppend1 = seismic(aa).headers(X1,:)./1000; 
    all_offset1 = cat(2,all_offset1,offsetToAppend1);
end  
%% this matrix includes offsets of all the required stations by removing first repetitive station 
all_offset_final = all_offset1(:,length(off1)+1:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%
allData = seismic.traces;
for i = 1:Number_of_segy_files
    dataToAppend = seismic(i).traces ;
    allData = cat(2,allData,dataToAppend);
end  
  
alldata_final = allData(:,length(off1)+1:end);
zero_matrix = find(all_offset_final==0);
alldata_final(:,zero_matrix)=0;
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
%% AUTO-CORRELATION       
[Dauto,t_cross]=auto_correlation_map(alldata_final,max_lag,dt);
Dauto(:,1:length(all_offset_final)) =  Dauto(:,1:length(all_offset_final))/SCALE_FACTOR_first_instrument;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% wiggle dipaly of auto-corelations%%%%%%%%%%%%%%%%%%%%%%
maximums = max(Dauto,[],'all');
figure 
pimage(all_offset_final,t_cross,Dauto)  
caxis([-maximums/contrast_coff maximums/contrast_coff]);
x=10;y=10;%starting screen position
     w=450;%figure width
     h=300;%figure hieght
    set(gcf,'position',[x y w h]);
    set(gca,'TickDir','out');
    set(gca,'xaxislocation','top');
    set(gca,'TickDir','out');
    ax = gca;
    ax.XAxis.LineWidth = 7;
    ax.YAxis.LineWidth = 7;
    fs=14;
%     xticks(-25:5:25);
%     yticks(0:0.1:1);
     xlim([4 4.3])
%     ylim([0 1])
    set(gca,'fontsize', fs,'FontWeight','Bold')
    set(gca,'FontName','Times New Roman')
    xlabel('Offset [km]', 'FontName','Times New Roman','FontSize', fs)
    ylabel('Time [s]', 'FontName','Times New Roman','FontSize', fs)
    lab=strcat('Auto-correlation');
    set(gcf, 'PaperPositionMode', 'auto' );
    print(gcf,'-r300','-dpng',lab)
    title('AUTO-CORRELATION')
%%%%%%%%%%%%%%%%%%%%%%%%
%% amplitude VS time for all the traces of all the OBS's
figure 
plot(t_cross,Dauto)
x=10;y=10;%starting screen position
     w=450;%figure width
     h=220;%figure hieght
    set(gcf,'position',[x y w h]);
    set(gca,'TickDir','out');
    set(gca,'xaxislocation','top');
    set(gca,'TickDir','out');
    ax = gca;
    ax.XAxis.LineWidth = 7;
    ax.YAxis.LineWidth = 7;
    fs=14;
%     xticks(-25:5:25);
%     yticks(0:0.1:0.5);
%     xlim([-25 25])
%     ylim([0 0.5])
    set(gca,'fontsize', fs,'FontWeight','Bold')
    set(gca,'FontName','Times New Roman')
    xlabel('Time [s]', 'FontName','Times New Roman','FontSize', fs)
    ylabel('Amplitude', 'FontName','Times New Roman','FontSize', fs)
    lab=strcat('amplitude spectrum_all_traces');
    set(gcf, 'PaperPositionMode', 'auto' );
    print(gcf,'-r300','-dpng',lab)
    title('Amplitude vs Time')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% amplitude VS time for summation of all the traces (source wavelet for positivie lag time)
DDauto=sum(Dauto,2);
figure
plot(t_cross,DDauto,'Color','r','LineWidth',2)
x=10;y=10;%starting screen position
     w=450;%figure width
     h=220;%figure hieght
    set(gcf,'position',[x y w h]);
    set(gca,'TickDir','out');
    set(gca,'xaxislocation','top');
    set(gca,'TickDir','out');
    ax = gca;
    ax.XAxis.LineWidth = 7;
    ax.YAxis.LineWidth = 7;
    fs=14;
%     xticks(-25:5:25);
%     yticks(0:0.1:0.5);
%     xlim([-25 25])
%     ylim([0 0.5])
    set(gca,'fontsize', fs,'FontWeight','Bold')
    set(gca,'FontName','Times New Roman')
    xlabel('Time [s]', 'FontName','Times New Roman','FontSize', fs)
    ylabel('Amplitude', 'FontName','Times New Roman','FontSize', fs)
    lab=strcat('Source wavelet');
    set(gcf, 'PaperPositionMode', 'auto' );
    print(gcf,'-r300','-dpng',lab)
    title('Source wavelet(Amplitude vs Time after summation)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% display of zero phase wavelet in terms of negative and positive lag time
%% (source wavelet for positive and negative lag time)
figure
DDauto_symmetric_neg = flip(DDauto(1:length(DDauto)/2));
DDauto_symmetric_pos = DDauto(1:length(DDauto)/2-1);
DDauto_symmetric = vertcat(DDauto_symmetric_neg,DDauto_symmetric_pos);
t_new = (-length(DDauto_symmetric)/2:(length(DDauto_symmetric)/2)-1)*dt;
plot(t_new,DDauto_symmetric,'Color','r','LineWidth',2)
x=10;y=10;%starting screen position
     w=450;%figure width
     h=220;%figure hieght
    set(gcf,'position',[x y w h]);
    set(gca,'TickDir','out');
    set(gca,'xaxislocation','top');
    set(gca,'TickDir','out');
    ax = gca;
    ax.XAxis.LineWidth = 7;
    ax.YAxis.LineWidth = 7;
    fs=14;
%     xticks(-25:5:25);
%     yticks(0:0.1:0.5);
%     xlim([-25 25])
%     ylim([0 0.5])
    set(gca,'fontsize', fs,'FontWeight','Bold')
    set(gca,'FontName','Times New Roman')
    xlabel('Time [s]', 'FontName','Times New Roman','FontSize', fs)
    ylabel('Amplitude', 'FontName','Times New Roman','FontSize', fs)
    lab=strcat('Source wavelet (positive-negative time)');
    set(gcf, 'PaperPositionMode', 'auto' );
    print(gcf,'-r300','-dpng',lab)
    title('Source wavelet (positive-negative time)')
    %%%%%%%%%%%%%%%%%%%%%%
    max_source = max(DDauto);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% exporting source file as output 
Source_matrix = DDauto;
zero_matrix= zeros(length(t)-length(DDauto),1);
Source_file = cat(1,Source_matrix,zero_matrix); 
Source_file((max_sample_time+1:end),:) = [];   
destination=strcat('/Volumes/DRIVE_DATA/TEST3/SOURCE_OUTPUT_FILE/eswir_shot_file');
dlmwrite(destination, Source_file, 'precision', 6, 'delimiter', '\t')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% power spectral density (PSD)
nt = length(t);
f=1/dt;
N = length(DDauto);
xdft = fft(DDauto);
xdft = xdft(1:N/2+1);
psdx = (1/(f*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
psdx=psdx/max(psdx);
freq = 0:f/length(DDauto):f/2;
figure,
plot(freq,(psdx),'b','LineWidth',6)
grid on
title('Power Spectral Density (PSD)')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
grid
x=10;y=10;%starting screen position
     w=450;%figure width
     h=220;%figure hieght
    set(gcf,'position',[x y w h]);
    set(gca,'TickDir','out');
    set(gca,'xaxislocation','top');
    set(gca,'TickDir','out');
    ax = gca;
    ax.XAxis.LineWidth = 7;
    ax.YAxis.LineWidth = 7;
    fs=14;
%     xticks(-25:5:25);
%     yticks(0:0.1:0.5);
%     xlim([-25 25])
%     ylim([0 0.5])
    set(gca,'fontsize', fs,'FontWeight','Bold')
    set(gca,'FontName','Times New Roman')
%     xticks(0:5:f);
%     xlim([0 40])
    lab=strcat('power spectral density');
    set(gcf, 'PaperPositionMode', 'auto' );
    print(gcf,'-r300','-dpng',lab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% amplitude spectrum
plot_spectral_attributes_amp(t_new(1,1),DDauto_symmetric,dt,f/2,1);
grid on
x=10;y=10;%starting screen position
     w=450;%figure width
     h=220;%figure hieght
    set(gcf,'position',[x y w h]);
    set(gca,'TickDir','out');
    set(gca,'xaxislocation','top');
    set(gca,'TickDir','out');
    ax = gca;
    ax.XAxis.LineWidth = 7;
    ax.YAxis.LineWidth = 7;
    fs=14;
%     xticks(-25:5:25);
%     yticks(0:0.1:0.5);
%     xlim([-25 25])
%     ylim([0 0.5])
    set(gca,'fontsize', fs,'FontWeight','Bold')
    set(gca,'FontName','Times New Roman')
%     xticks(0:2.5:f);
%     xlim([0 20])
    lab=strcat('amplitude spectrum');
    set(gcf, 'PaperPositionMode', 'auto' );
    print(gcf,'-r300','-dpng',lab)
    title('Amplitude Spectrum')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Important: for display of phase spectrum, the symmetric matrix of autocoorelation must be used
plot_spectral_attributes_phase(t_new(1,1),DDauto_symmetric,dt,f/2,1);
   
grid on
x=10;y=10;%starting screen position
     w=450;%figure width
     h=220;%figure hieght
    set(gcf,'position',[x y w h]);
    set(gca,'TickDir','out');
    set(gca,'xaxislocation','top');
    set(gca,'TickDir','out');
    ax = gca;
    ax.XAxis.LineWidth = 7;
    ax.YAxis.LineWidth = 7;
    fs=14;
    set(gca,'fontsize', fs,'FontWeight','Bold')
    set(gca,'FontName','Times New Roman')
%     xticks(0:2.5:f);
%     xlim([0 20])
%     yticks(-180:90:180);
    lab=strcat('phase spectrum');
    set(gcf, 'PaperPositionMode', 'auto' );
    print(gcf,'-r300','-dpng',lab)
    title('Phase Spectrum')
