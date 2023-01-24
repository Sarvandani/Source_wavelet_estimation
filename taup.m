function data_out=taup(data_in,x,p,dt,nt,wn,ind,verbose)
%
%     function data_out=taup(data_in,x,p,dt,nt,wn,ind,verbose)
%
%     A forward and inverse tau-p transform function.
%     data_in is the input array tau/time down columns, x/p along rows
%     x is the vector of receiver locations. Note dx assumed constant.
%     p is the vector of slownesses
%     dt is the time series sample rate
%     nt is the number of time or tau samples
%     wn is the faction of white noise to stabilize matrix e.g. 0.02
%     data_out is the output transformed array
%     ind=+1: forward tau-p from x-t to p-tau domain
%     ind=-1: inverse tau-p from p-tau to x-t domain
%     verbose=0; no status messages, 1 status messages
%
%     Written by Richard T. Coates 14 Feb 95 using eqns (7) and (8) from
%     Foster & Mosher, 1992. Suppression of multiple reflections using 
%     radon transform Geophys, 57, 386-395.
%

%     set up constants
  nt=2*floor((nt+0.1)/2);
  ci=0+i;
  dw=2*pi./(nt*dt);

% orientate x and p vectors
    [n,m]=size(x);
    if( n > m );x=x';end;
    nx=max(n,m);
    [n,m]=size(p);
    if( n < m );p=p';end;
    np=max(n,m);
    dx=x(2)-x(1);
    
% forward transform from space time to slowness tau
  if(ind == 1)

% do fft over time
    data_in=fft(data_in,nt);
    data_out=zeros(np,nt);
    
% for each frequency build y=R.data, RR* and invert
    for j=1:nt/2+1;

      if(verbose == 1);
	if(mod(j,50) == 0);
	  disp('Doing frequency number');disp(j);
	end;
      end;

      r=dx*exp(ci*(j-1)*dw*p*x);
      a=r*r';
      w=wn*max(abs(diag(a)));
      a=a+w*diag(ones(1,np));
      y=r*conj(data_in(j,:)');
      data_out(:,j)=a\y;

    end;    
    for j=2+nt/2:nt
      data_out(:,j)=conj(data_out(:,nt+2-j));
    end;
    
% do inverse fourier transform
    for j=1:np;
       data_out(j,:)=ifft(data_out(j,:));
    end;

% transpose real results
    data_out=real(data_out');

  else
    
% do fft over time
    data_in=fft(data_in,nt);
    data_out=zeros(nx,nt);

% for each frequency build R and form y=R.data,     
    for j=1:nt/2+1

      r=dx*exp(ci*(j-1)*dw*p*x);
      data_out(:,j)=r'*conj(data_in(j,:)');      
      if(verbose == 1);
	if(mod(j,50) == 0);
	  disp('Doing frequency number');disp(j);
	end;
      end;


    end;    
    for j=2+nt/2:nt
      data_out(:,j)=conj(data_out(:,nt+2-j));
    end;
    
% do inverse fourier transform
    for j=1:nx;
       data_out(j,:)=ifft(data_out(j,:));
    end;

% transpose real results
    data_out=real(data_out)';
    
  end;    
