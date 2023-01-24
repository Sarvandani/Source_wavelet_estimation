function impedance=s_impedance(refl,scale)
% Function computes impedance from reflection coefficients
% Start time and sample interval of the output impedance is that
% of the reflection coefficients
% Method used: Recursive computation from reflectivity
% See also: s_refl2impedance
%
% Written by: E. R.: February 14, 2002
% Last update: April 1, 2002: Handle multi-trace input
%
%           impedance=s_impedance(refl,scale)
% INPUT 
% refl      seismic dataset with reflection coefficients
% scale     Scale information. Either a seismic structure with an impedance
%           trace or a scalar representing the impedance at the first sample
%           of the output. In the former case the number of samples of 
%           scale.trace must equal 1 + the number of samples of refl.trace.   
% 
% OUTPUT
% impedance   seismic data set with computed impedance

impedance=refl;

ntr=size(refl.traces,2);
temp=(1+refl.traces)./(1-refl.traces);
temp=[ones(1,ntr);cumprod(temp)];
if isnumeric(scale)
   impedance.traces=scale*temp;
   htext=['First impedance value is ',num2str(scale)];

elseif isstruct(scale)
   if length(scale.traces) ~= length(temp)
      error(' Traces in "scale" and "impedance" do not have the same number of samples')
   end
   temp1=log(scale.traces./temp);
   idx=~isnan(temp1);
   a=exp(median(temp1(idx)));
   impedance.traces=a*temp;
   htext='Impedance scaled to reference';
else
   error(' Illegal variable type for "scale"')
end

%	Check for NaNs
if any(isnan(impedance.traces(:)))
   impedance.null=NaN;
end
impedance.last=refl.last+refl.step;

%    Append history field
if isfield(refl,'history')
   impedance=s_history(impedance,'append',htext);
end
