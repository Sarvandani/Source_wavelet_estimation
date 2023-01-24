function title_handle=mytitle(string,varargin)
% Function creates a "standard" title
% Written by: E. Rietsch: April 5, 2001
% Last updated: September 20, 2001: Handle underscores
%
%           mytitle(string,varargin)
% INPUT
% string    title string
% % varargin  one or more cell arrays; the first element of each cell array is a keyword,
%           the other elements are parameters. Presently, keywords are:
%      'fig'      figure handle. Default: {'fig',gcf} i.e. current figure
%      'color'    Title color. Default {'color','r'}
%      'fontsize' Font size of title text. Default: {'fontsize',ftsze}
%                 where ftsze=15*min([1,66/(length_of_title_string)]) for landscape format
%                       ftsze=15*min([1,48/(length_of_title_string)]) for portrait format
%      'fontname' Font name. Default: {'fontname','Arial'}

global S4M

% 	Set default parameters
param.fig=gcf;
param.color='r';
param.fontsize=[];
param.fontname=S4M.font_name;
param.tex='no';

%	Read/interpret input arguments
param=assign_input(param,varargin);

if isempty(param.fontsize)
  orient=get(param.fig,'PaperOrientation');
  if iscell(string)
    lstring=0;
    for ii=1:length(string)
      lstring=max([lstring,length(string{ii})]);
    end
    sph=gca;
    pos=get(sph,'Position');
    pos(2)=0.95*pos(2);
    pos(4)=0.91*pos(4);
    set(sph,'Position',pos); 
  else
    lstring=length(string);
  end
  if strcmp(orient,'landscape')
    param.fontsize=min([1,66/lstring])*15;
  else
    param.fontsize=min([1,40/lstring])*15;
  end
end

figure(param.fig);

if strcmpi(param.tex,'no')
   string=mnem2tex(string);
end

htitle=title(string,'Color',param.color,'FontSize',param.fontsize);

if isempty(param.fontname)
   set(htitle,'FontName',param.fontname)
end

if nargout == 0
   clear title_handle
end