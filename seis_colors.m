function RGB = seis_colors(n)
% This code provides a seismic color map for various seismic data gathers.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0 
    n = size(get(gcf,'colormap'),1);
end
m    = ceil(n/3);
Top  = ones(m,1);
Bot  = zeros(m,1);
up   = (0:m-1)'/m;
down = flipud(up);
R    = [Bot; up; 1; Top; down];
G    = [Bot; up; 1; down; Bot];
B    = [up; Top; 1; down; Bot];
RGB  = [R G B];
xlarge = 4*m+1-n;
xblue  = round(xlarge/2);
xred   = xlarge - xblue;
RGB([1:xblue 4*m-xred+2:4*m+1],:) = [];
