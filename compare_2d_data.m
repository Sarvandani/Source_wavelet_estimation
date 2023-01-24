% This function is to process and compare the observed and modeled data
path = '../../swir_2D_V4_synt601/';

dlen=0.03;

ishot=1; plottrace=[240:10:320];
ishot=3; plottrace=[10:20:140]+10;
ishot=2; plottrace=[280:10:320];

% near offset data
ishot=3; plottrace=[170:22:220 270:25:330];
ishot=2; plottrace=[40:30:120 180:20:222];
ishot=1; plottrace=[10:20:60 135:20:192];

ites=[0 1 25]

maxwig=160*3;
% maxwig=600*1;

dataplot=[0 1 0 0];
% dataplot=[0 0 0 0];
% options to plot data panel, [overlying wiggle; seperate wiggle; 
% residuals wiggle (1) or image (2), horizontally plot of selected traces]

ifspec=1;  % if call compare_data_function for more analysis
remute=1;
ifnorm=0;
ifreduce=1; 

ifformod=0;
ifrawdata=0; 

if_corr_airgun=0; % if correct the amplitude error due to the air-gun array problem
    % Do not use this in real6 series, sicne I have done it when convering
    % the data.
ifsecondlabel=0;  % the default xaxis is offset, the second choise is location 

if ifnorm>0 && ifspec>0 && remute==0
    disp('You are going to analyse the normalised data, are you sure?');
    keyboard;
end
if ifnorm==1; maxwig=15; end


if ~isempty(strfind(path,'synt'))
    prjn='synt';
elseif ~isempty(strfind(path,'kern'))
    prjn='kern';
else
    prjn='swir';
end


colorseq=['b' 'r' 'k' 'c' 'm'];
caxisv=[-1 1]*200;

% parameters for time axe
tsteps = 6001;
% tsteps = 4001;
dt = 0.002;

if path(strfind(path,'real')+4)>='5', dlen=0.03; 
elseif path(strfind(path,'real')+4)=='4', dlen=0.04; end

% dlen=0.015;

rec_time = (0:tsteps-1)' * dt*ones(1,length(ites));

% rcvr's specs
nrecs = 331;
% nrecs = 287;

filename_in = [path, './data/swir_V4_rcvr_2D'];
filename_in = [path, './data/swir_V4_rcvr_2Dreal'];
rcvr_in=load(filename_in,'-ascii');
filename_in = [path, './data/swir_V4_shot_2D'];
filename_in = [path, './data/swir_V4_shot_2Dreal_30m'];
% filename_in = [path, './data/swir_V4_shot_kernel40'];
shot_in=load(filename_in,'-ascii');
% xplotpos=(rcvr_in(:,2)-1)*dlen;
% xplotpos=1:nrecs;
xplotpos=(rcvr_in(:,2)-shot_in(ishot,3))*dlen;
xplotpos2=(rcvr_in(:,2)-1)*dlen;
xplotpos3=rcvr_in(:,1);
xplotpos3=xplotpos;

filename_in = sprintf('%s/data/%s_mut_%04dhydr',path,prjn,ishot);
mute_in=load(filename_in,'-ascii');
tw=sqrt(((shot_in(ishot,4)-rcvr_in(1,3))*0.001)^2+xplotpos.^2)/1.5...
    -0.75;
tw1=round(tw/dt);
topmute=mute_in(:,2); topmute=round(topmute);
% topmute(topmute<tw1-10)=topmute(topmute<tw1-10)-100;
botmute=mute_in(:,3); 
botmute=min([topmute+250 botmute],[],2);
% botmute=min([topmute+250 tw1],[],2);
% botmute=max([topmute+300 round(tw/dt)],[],2);

botmute((botmute-topmute)<200)=1;

% automatically generate the xrange and yrange for plotting data
itmp=find(topmute<botmute);
iltmp=min(itmp); irtmp=max(itmp);
xrange=xplotpos([iltmp irtmp])+...
    [-0.2 0.5]'*0.005*abs(diff(xplotpos([1 end])));
yrange=rec_time([max(1,min(topmute(itmp))) max(botmute(itmp))],1)+...
    [-1 1]'*0.2;

topmute_t=topmute*dt; 
botmute_t=botmute*dt; 
% botmute_t=nan;

if ifrawdata==1, 
    ites=[0]; yrange=[2 11];
    xrange=[min(xplotpos) max(xplotpos)];
    remute=0; ifreduce=0;
    dataplot=[0 2 0 0];
end

if ifreduce==1
    yrange=[1.2 3.3];
    vred=6.5;
end

% input data
nite=length(ites);
dmod=nan(tsteps,nrecs,nite);
dmod_mute=dmod;

if dataplot(1)>0, fig_overlay=figure('position',[2900 500 850 360*nite],...
        'defaultaxesfontname','helvetica','defaultaxesfontsize',18); end
if dataplot(2)>0,fig_wig=figure('position',[1900 500 650 360*nite],...
        'defaultaxesfontname','helvetica','defaultaxesfontsize',18); end
dratio=[2.5 1 1];
dratio=[0];

for i=1:nite
    ite=ites(i);
    if ifformod==1
        if ite==0
            fdata = sprintf('%s/data/%s_ffid_%04dhydr',path,prjn,ishot);
%             fdata = sprintf('%s/data/synt_ffid_%04dhydr',path,ishot);
        else
%             fdata = sprintf('%s/dobs.mplx.hydr.s%03d.ite%03d',path,ishot,ite);
            fdata = sprintf('%s/dobs.mplx.hydr.s%03d',path,ishot);
        end
    elseif ifrawdata==1
        if ite==0
            fdata = sprintf('%s/data/%s_ffid_%04dhydr',path,prjn,ishot);
        end
    else
        if ite==0
            fdata = sprintf('%s/dobs.mplx.hydr.s%03d',path,ishot);
            fdata = sprintf('%s/data/%s_ffid_%04dhydr',path,prjn,ishot);
        else
            fdata = sprintf('%s/dmod.hydr.s%03d.ite%03d',path,ishot,ite);
        end
    end
    fdata
    fid = fopen (fdata, 'r');
    [data_in] = fread (fid, Inf, 'float');
    fclose(fid);
    
    dmod(:,:,i)=reshape(data_in,[nrecs tsteps])';
    
%     if i==1;  dmod(:,:,i)=dmod(:,:,i)/1; end
    
    if remute==1
        dmod_mute(:,:,i)=mute_data(squeeze(dmod(:,:,i)),topmute,...
            botmute);
        dmod(:,:,i)=dmod_mute(:,:,i);
    end
    
    if ifnorm==1
        dmod(:,:,i)=normalise_data_amplitude(dmod(:,:,i),0);
    end
    
    if ifreduce==1
        dmod(:,:,i)=reduce_data(dmod(:,:,i),xplotpos,vred,dt);
    end
    
    if dataplot(1)>0
        figure(fig_overlay);
        subplot(nite,1,i);
        if i==1
            seisplot_wigb_pcolor(dataplot(1),dmod(:,:,[i nite]),...
                xplotpos,rec_time(:,1),maxwig);
            set(gca,'xlim',xrange,'ylim',yrange);
            if dratio(1)>0, daspect(dratio); end
            if ifreduce==1
%                 ylabel(sprintf('Reduced time (s) at %3.1fkm/s',vred),...
%                     'Fontsize',18);
                ylabel('Reduced time (s) ');
            end
            %             legend(num2str([1 nite]));
        else
            seisplot_wigb_pcolor(dataplot(1),dmod(:,:,[i 1]),...
                xplotpos,rec_time(:,1),maxwig);
%             seisplot_wigb_pcolor(dataplot(1),dmod(:,:,[1 i]),...
%                 xplotpos,rec_time(:,1),maxwig);
            set(gca,'xlim',xrange,'ylim',yrange);
            if dratio(1)>0, daspect(dratio); end
            if ifreduce==1
%                 ylabel(sprintf('Reduced time (s) at %3.1fkm/s',vred),...
%                     'Fontsize',18);
                ylabel('Reduced time (s) ');
            end
%             legend(num2str([i 1]));
        end
        title(sprintf('Shot %d: Input dataset, ite %d  ',ishot,ite));
    end
    if dataplot(2)>0
        figure(fig_wig);
        subplot(nite,1,i);
        seisplot_wigb_pcolor(dataplot(2),dmod(:,:,i),...
            xplotpos,rec_time(:,1),maxwig);
        set(gca,'xlim',xrange,'ylim',yrange);
        if dratio(1)>0, daspect(dratio); end
        if ifreduce==1
%             ylabel(sprintf('Reduced time (s) at %3.1fkm/s',vred),...
%                 'Fontsize',18);
            ylabel('Reduced time (s) ');
        end
        plot(xplotpos,topmute_t-abs(xplotpos)/vred,'r','linewidth',2);
%         plot(xplotpos,topmute_t-abs(xplotpos)/vred+0.6,'b','linewidth',2);
        plot(xplotpos,botmute_t-abs(xplotpos)/vred,'b','linewidth',2);
        set(gca,'xtick',[-60:2:60]);
        xlabel('Source-receiver offset (km) ');
        if ifsecondlabel==1
            if ishot==1
                set(gca,'xtick',[-36:2:36]+18,'xticklabel',num2str([-36:2:36]'));
            elseif ishot==2
                set(gca,'xtick',[-36:2:36]+10.2,'xticklabel',num2str([-36:2:36]'));
            elseif ishot==3
                set(gca,'xtick',[-36:2:36]-13.8,'xticklabel',num2str([-36:2:36]'));
            end
            xlabel('Across-axis distance (km) ');
        end
        title(sprintf('Input dataset, ite %d  ',ite));
    end
end

if dataplot(3)>0
    figure('position',[1900 500 850 400*nite],...
        'defaultaxesfontname','helvetica','defaultaxesfontsize',18);
    for i=1:nite
        ite=ites(i);
        subplot(nite,1,i);
        if i==1            
            seisplot_wigb_pcolor(dataplot(3),dmod(:,:,i),...
                xplotpos,rec_time(:,1),maxwig);
            set(gca,'xlim',xrange,'ylim',yrange);
            plot(xplotpos,topmute_t-abs(xplotpos)/vred,'r','linewidth',2);
            plot(xplotpos,botmute_t-abs(xplotpos)/vred,'b','linewidth',2);
            if dratio(1)>0, daspect(dratio); end
            if ifreduce==1
%                 ylabel(sprintf('Reduced time (s) at %3.1fkm/s',vred),...
%                     'Fontsize',18);
                ylabel('Reduced time (s) ');
            end
            title(sprintf('Shot %d:Observed data, ite %d  ',ishot,ite));
        else
            seisplot_wigb_pcolor(dataplot(3),diff(dmod(:,:,[i,1]),[],3),...
                xplotpos,rec_time(:,1),maxwig);
            if ifreduce==1
%                 ylabel(sprintf('Reduced time (s) at %3.1fkm/s',vred),...
%                     'Fontsize',18);
                ylabel('Reduced time (s) ');
            end
            set(gca,'xlim',xrange,'ylim',yrange);
            plot(xplotpos,topmute_t-abs(xplotpos)/vred,'r','linewidth',2);
            plot(xplotpos,botmute_t-abs(xplotpos)/vred,'b','linewidth',2);
            if ifsecondlabel==1
                if ishot==1
                    set(gca,'xtick',[-36:2:36]+18,'xticklabel',num2str([-36:2:36]'));
                elseif ishot==2
                    set(gca,'xtick',[-36:2:36]+10.2,'xticklabel',num2str([-36:2:36]'));
                elseif ishot==3
                    set(gca,'xtick',[-36:2:36]-13.8,'xticklabel',num2str([-36:2:36]'));
                end
                xlabel('Across-axis distance ');
            end
            if dratio(1)>0, daspect(dratio); end
            title(sprintf('Difference to observed data, ite %d  ',ite));
        end
    end
end

% scalefac=5;
% if ifformod==1
% %     dmod(:,:,1)=dmod(:,:,1)*scalefac;
%     dmod(:,:,2)=dmod(:,:,2)/scalefac;
%     if remute==1
% %         dmod_mute(:,:,1)=dmod_mute(:,:,1)*scalefac;
%         dmod_mute(:,:,2)=dmod_mute(:,:,2)/scalefac;
%     end
% end

% Plot
fontsize=18;
titlesize=20;

if dataplot(4)==1
    nplot=length(plottrace);
    traces=nan(tsteps,nite); traces_mute=traces;
    ampfac=nan(nplot);
    fig1=figure('Position',[1500 100 1200 400*nplot]);
    % fig2=figure('Position',[2200 100 1200 400*nplot]);
    for iplot=1:nplot
        itra=plottrace(iplot);
        traces=squeeze(dmod(:,itra,:));
        traces_mute=squeeze(dmod_mute(:,itra,:));
        stdobs=std(traces_mute(:,1));
        stdmod=std(traces_mute(:,2));
        
        %     stdobs=std(traces(:,1));
        %     stdmod=std(traces(:,2));
        %     ampfac(iplot)=stdobs/stdmod;
        ampfac(iplot)=1;
        traces(:,2:end)=traces(:,2:end)*ampfac(iplot);
        traces_mute(:,2:end)=traces_mute(:,2:end)*ampfac(iplot);
        
        figure(fig1);
        subplot(nplot,1,iplot);
        if remute==1
            plot(rec_time,traces_mute);
        else
            plot(rec_time,traces);
        end
        legend(num2str(ites'));
        grid on; hold on;
        title(sprintf('Trace %3d, amp factor %7.6f ',itra,ampfac(iplot)));
        %     set(gca,'ylim',[-1 1]*max(abs(traces(:))));
        
        %     figure(fig2);
        %     subplot(nplot,1,iplot);
        %     plot(rec_time,diff(traces,1,2));
        %     legend(num2str(ites'));
        %     grid on; hold on;
        %     set(gca,'fontsize',fontsize);
        %     title(sprintf('Trace %3d, amp factor %7.6f',itra,ampfac(iplot)),...
        %         'Fontsize',titlesize);
        %     set(gca,'xlim',...
        %         [-3 0]+min(mute(itra,3)*dt+1,10),...
        %         'ylim',[-1 1]*max(abs(traces_mute(:))));
        
        %     keyboard
    end
end

if ifspec==1
    if remute
        datain=dmod_mute(:,:,:);
    else
        datain=dmod(:,:,:);
    end
    
    if  if_corr_airgun==1
        fprintf('The air-gun errors will be corrected for shots ID > 270\n');
        airgun_corr_amp_fac=interp1([0 270 280 500],[1 1 4/3 4/3],rcvr_in(:,1));
        for irec=1:nrecs
            datain(:,irec,1)=datain(:,irec,1)*airgun_corr_amp_fac(irec);
        end
    end
    compare_data_function
end
