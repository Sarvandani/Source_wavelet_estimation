% %% Define mute file
clear all
clc
max_num_shot = 182
%I made this cofficent in the code. intitial value was 35,
%we should change it based on the data. if this value is larger
% we will have less amount of data
mute_coffic = 35
%another important thing is to check seusmic.headers.
X = 4
%% Define mute file
for nn=1:1
    lab=strcat('0',num2str(nn),'b_obs',num2str(nn),'.segy');
    SEISMIC=read_segy_file(lab);
    data=SEISMIC.traces(:,1:length(SEISMIC.traces(1,:)));
    data_cut=data;
    mut_top=[];
    mut_bot=[];
    for I=1:length(data(1,:))
        %find with ~= shows the number of 1s replaced with the numbers were not equal to 0.1 value.
        J=find(data(:,I)~=0,1);
        if isempty(J)~=1
        %any number deducted from j will shift the top mute on shot gather
            mut_top(I)= J+50;
        else
            mut_top(I)=mut_top(I-1);
        end
        header(:,:) = SEISMIC.headers;
        ind=round(abs(SEISMIC.headers(X,(I)))/mute_coffic);
        mut_bot(I)=mut_top(I)+ind;
        %if we remove this part, the mute line and area are not compatible
        data_cut(mut_bot(I):end,I)=0;
        data_cut(1:mut_top(I),I)=0;
    end
    %find with >= %find with ~= shows the number of 1s in a matrix replaced with the numbers were bigger than 0.1 value.
    %nnn=find(SEISMIC.headers(4,:)>=0,1);
    %mut_bot(nnn-15:nnn+20)=mut_top(nnn-15:nnn+20)-95;
    close all
    imagesc(data);
    caxis([-1e5 1e5])
    colormap gray
    hold on
    plot(1:max_num_shot,mut_top,'b')
    plot(1:max_num_shot,mut_bot,'b')
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
    
    %off=SEISMIC.headers(4,:)./1000;
    %t=0:(SEISMIC.step)./1000:(SEISMIC.last)./1000;
    % this close all prevents from overlapping the figures
    close all
    %fs=16;
    % I don't think if we need this part
    for mm=1:length(SEISMIC.traces(1,:))
        SEISMIC.traces(1:mut_top(mm)*2,mm)=0;
        SEISMIC.traces(mut_bot(mm)*2:end,mm)=0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Yas:This part is only for plotting the gathers
    %Yas:we may need this part for frq filtering
    %SEISMIC=s_filter(SEISMIC,{'ormsby',[1 3 10 25]});
    %figure
     hold on
     set(0,'defaultfigurecolor',[1 1 1])
    off=SEISMIC.headers(X,:)./1000;
    t=(0:SEISMIC.step/1000:SEISMIC.last/1000);
    %Yas:mwigb is used for trace display on the figure, without it we 
    % don't have any figure.
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
    % Yas:this part is used to play only mute lines
    plot(-1*(off),SEISMIC.step/1000* 2*(mut_top),'r')
    plot(-1*(off),SEISMIC.step/1000* 2*(mut_bot),'r')
    set(gca,'fontsize', fs)
    set(gca,'FontName','Times New Roman')
    xlabel('Signed Offset [km]', 'FontName','Times New Roman','FontSize', fs)
    ylabel('Time [s]', 'FontName','Times New Roman','FontSize', fs)
    lab=strcat('Fig obs', num2str(nn),'_input_with_mutes');
    set(gcf, 'PaperPositionMode', 'auto' );
    print(gcf,'-r300','-dpng',lab)
    
%     KK=find(all(ismember(SEISMIC.header_info,'ffid'),4)==1);
%     SEISMIC.headers(KK,:)=nn;
%     lab1=strcat('eswir_ffid_',num2str(nn),'hydr.segy');
%     write_segy_file(SEISMIC,lab1);
end
