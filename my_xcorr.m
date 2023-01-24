function rxx=my_xcorr(x,N,d)
% This code is written to x-correlation of a 1xN seismic trace x with a vector d. 
% The required inputs are:
% x: seismic shot gather(s)
% N: number of lags
% d: vector (if not provided then the function correlates the vector x)
% 
% The output is:
% rxx: the x-correlation vector output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code if for the book titled: Processing Surface Seismic Reflection
% Data Using MATLAB by Wail A. Mousa & Abdullatif A. Al-Shuhail
% September 2011.

if nargin==2
    rxx=zeros(1,N);
    L=length(x);
    for n=1:N
        for l=1:L
            if l+n <= L+1
                rxx(n)=rxx(n)+x(l)*x(l+n-1);
            end
        end
    end

else

    rxx=zeros(1,N);
    L=max(length(d),length(x));
    if length(d)>length(x)
        x=[x,zeros(1,L-length(x))];
    elseif length(d)<length(x)
        d=[d,zeros(1,L-length(d))];
    end

    for n=1:N
        for l=1:L
            if l+n <= L+1
                rxx(n)=rxx(n)+x(l)*d(l+n-1);
            end
        end
    end
end
    
