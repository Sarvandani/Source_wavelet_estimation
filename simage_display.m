function simage_display(D,x,z,s)
% This code is written to display gray and colored variable density displays 
% The required inputs are:
% D: seismic data
% x: Offset or trace number values
% z: time or depth axis values
% s: must be 0 for gray and must be 1 for color displays 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code if for the book titled: Processing Surface Seismic Reflection
% Data Using MATLAB by Wail A. Mousa & Abdullatif A. Al-Shuhail
% September 2011.

if s==0 %Gray scale variable density display
    pcolor(x,z,-1*D),
    shading interp;
    axis ij;
    colormap(gray);%Using a gray-scaled colormap
    colorbar
    % Displaying the seismic gather in the traditional way for displaying seismic data
    x=10;y=10;%starting screen position
    w=300;%figure width
    h=600;%figure hieght
    set(gcf,'position',[x y w h]);
    % Putting the spatial axis in the top of the figure 
    set(gca,'xaxislocation','top')
elseif s==1
    pcolor(x,z,-1*D),
    shading interp;
    axis ij;
    colormap(seis_colors);%Using a colored colormap
    colorbar
    % Displaying the seismic gather in the traditional way for displaying seismic data
    x=10;y=10;%starting screen position
    w=300;%figure width
    h=600;%figure hieght
    set(gcf,'position',[x y w h]);
    % Putting the spatial axis in the top of the figure 
    set(gca,'xaxislocation','top')       
end