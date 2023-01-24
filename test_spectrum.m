clear all
clc
%% manual and automatic periodgrams are not the same when we have a matrix
%% for a vector they are the same.
X = 3
%% this chnage the range of amplitude and make different contrast on the figure
contrast_coff = 1

for nn=2:2
    lab=strcat('0',num2str(nn),'b_obs',num2str(nn),'.segy');
    seismic=read_segy_file(lab);
         seismic=s_resample(seismic,2);

    off=seismic.headers(X,:)./1000;
    t=(0:seismic.step/1000:seismic.last/1000);
    data=seismic.traces(:,1:length(seismic.headers(1,:)));
    %%%%%%%%%%%%%%%
%      rng default
% n = 0:999;
% data = cos(pi/4*n) + randn(size(n));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This method is for db but I dont know why my max freq is two times more
figure,
     s_spectrum(seismic,{'plot','amp'},{'timezero','actual'},{'figure','old'}, {'scale','db'},...,
     {'frequencies',0,250},{'normalize','no'}, {'colors','k'})

%     auto=randn(100,1);
%     s_spectrum(data,{'plot','amp'},{'timezero','actual'},{'figure','old'}, {'scale','db'},...,
%     {'frequencies',0,25},{'normalize','yes'}, {'colors','k'})
    
    dt = 0.002;
    nt = 30001;
 fs=1/dt;
 %% This method (Wail A. Mousa and Abdullatif A. Al-Shuhail 2011) was for averagine trace of seimsic data, and we should not take data1/max(data1) for other dataest
% [data1,f] = periodogram(data,rectwin(length(data)),length(data),fs);
% data1=data1/max(data1);
% data1=20*log10(abs(data1));
% 
% figure,plot(f,data1,'r--')
% xlabel('Frequency (Hz)','FontSize',14)
% ylabel('Normalized PSD','FontSize',14)
% grid
% legend('Before decon','After decon')
%%%%%%%%%%%%%%%%%%%%
%% manual periodgram db/hz
%% https://fr.mathworks.com/help/signal/ug/power-spectral-density-estimates-using-fft.html
N = length(data);
xdft = fft(data);
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
psdx=psdx/max(psdx);
freq = 0:fs/length(data):fs/2;
figure,
plot(freq,10*log10(psdx))
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% automatic periodgram (db/hz)
figure, 
periodogram(data,rectwin(length(data)),length(data),fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% amplitude in freq domain
%% 
plot_spectral_attributes(min(data),data,dt,250,1);
    
end
