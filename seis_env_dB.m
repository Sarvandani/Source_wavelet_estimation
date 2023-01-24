function seis_env_dB(Data,Data_g,t,tnum)
% This code is written to the amplitude envelope for seismic traces before & after gain control. The required inputs are:
% Data: seismic shot gather(s) before gain control
% Data_g: seismic shot gather(s) after gain control
% t: the time axis values
% tnum: trace number
% 
% The outputs are one amplitude envelope in dB figures: either for the
% selected trace or for the average traces.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code if for the book titled: Processing Surface Seismic Reflection
% Data Using MATLAB by Wail A. Mousa & Abdullatif A. Al-Shuhail
% September 2011.



if nargin==4
    trace1=Data(:,tnum)';
    trace1=trace1/max(trace1);
    trace1_env_dB=20*log10(abs(hilbert(trace1+eps)));%log of envelop
    trace1_towp=Data_g(:,10)';%topw
    trace1_towp=trace1_towp/max(trace1_towp);
    trace1_env_dB_towp=20*log10(abs(hilbert(trace1_towp+eps)));%log of envelop %check envelop
    figure
    plot(t,trace1_env_dB,'b',t,trace1_env_dB_towp,'r-.','LineWidth',3)
    xlabel('Time (s)','FontSize',14)
    ylabel(['Envelope (dB) trace number :',num2str(tnum),''],'FontSize',14)
    grid
    legend('Before gain correction','After gain correction',3)
else
    trace1=mean(Data');
    trace1=trace1/max(trace1);
    trace1_env_dB=20*log10(abs(hilbert(trace1+eps)));%log of envelop
    trace1_towp=mean(Data_g');%topw
    trace1_towp=trace1_towp/max(trace1_towp);
    trace1_env_dB_towp=20*log10(abs(hilbert(trace1_towp+eps)));%log of envelop %check envelop
    figure
    plot(t,trace1_env_dB,'b',t,trace1_env_dB_towp,'r-.','LineWidth',3)
    xlabel('Time (s)','FontSize',14)
    ylabel('Envelope (dB) of the average trace','FontSize',14)   
    grid
    legend('Before gain correction','After gain correction',3)
end

