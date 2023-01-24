function [v_stack,t_stack]=vel_picking(Dsort,Hsort,vmin,dv,nv,cmp_start,cmp_end,cmp_step,num_points)
% This code is written to perform nmo/stacking velocity picking from CMP seismic gathers' sembelances.
% The required inputs are:
% Dsort: the sorted data
% Hsort: the sorted header
% vmin: minimum expected velocity
% dv: velocity sampling interval
% nv: number of velocity values to be used
% cmp_start: First CMP gather 
% cmp_end: Last CMP gather 
% cmp_step: CMP gather step size to be used
% num_points: number of velocity points to be picked and stored
% 
% The output is:
% v_stack: the picked velocity values in a vector
% t_stack: the corresponding times in a vector of the same length as
% v_stack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code if for the book titled: Processing Surface Seismic Reflection
% Data Using MATLAB by Wail A. Mousa & Abdullatif A. Al-Shuhail
% September 2011.


vmax =(nv-1)*dv+vmin;
cmp_num=cmp_start:cmp_step:cmp_end;
v_stack=zeros(length(cmp_num),num_points);%Picking 3 points
t_stack=v_stack;
R = 2;
L = 10;

for kkk=1:length(cmp_num)
    disp(['I am CMP number: ',num2str(cmp_num(kkk)),''])
    [Da,dt,t,cdp,jj,cmp_offset1]=extracting_cmp(Dsort,Hsort,cmp_num(kkk));
    [nt,nx]=size(Da);
    h =cmp_offset1;
    [S,tau,v] = velan(Da,dt,h,vmin,vmax,nv,R,L);
    
    %Displaying the data in wiggle 
    scale=1;
    if nx>1;
        %Raw data
        figure,mwigb(Da,scale,cdp,t)
        xlabel('Trace number','FontSize',14)
        title(['CMP: ',num2str(cmp_num(kkk)),''],'FontSize',14)%Note this part
        ylabel('Time (s)','FontSize',14)
        
        
    else
        plot(Da,t,'k'),
        hold on
        area(Da.*(Dcmp>0),t),colormap([0,0,0]),
        hold off
        axis tight
        xlabel(['CMP gather: ',num2str(cmp_num(kkk)),''],'FontSize',14)%Note this part
        ylabel('Time (s)','FontSize',14)
        set(gca,'YDir','reverse','XAxisLocation','Top')
    end
        
    figure,
    pcolor(v,tau,S);shading interp;axis ij;
    grid 
    xlabel('V (ft/s)','FontSize',14);
    ylabel('Time (s)','FontSize',14);
    colormap(jet),colorbar
    x=10;y=10;%starting screen position
    w=300;%figure width
    h=600;%figure hieght
    set(gcf,'position',[x y w h]);
    % Putting the spatial axis in the top of the figure 
    set(gca,'xaxislocation','top')


    
    %Picking: Define t-v pais for NMO
    %[v_stack(kkk,:),t_stack(kkk,:)]=ginput(num_points);   
    hold on
    % Loop, picking up the points.
    for i = 1:num_points
        [v_stack(kkk,i),t_stack(kkk,i)] = ginput(1);
        plot(v_stack(kkk,i),t_stack(kkk,i),'m+','LineWidth',2)% ko is a circlur point. you can use '*' point.
    end
    
    
    

end