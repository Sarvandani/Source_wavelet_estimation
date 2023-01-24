function stacking_chart(sx,gx,num_shots,num_trace_per_sg)
% This code is written to display a stacking chart for the given seismic data parameters.
% The required inputs are:
% sx: sources x-axis locations
% gx: receivers x-axis locations
% num_shots: is the number of shot gathers
% num_trace_per_sg: number of traces/shot gather
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code if for the book titled: Processing Surface Seismic Reflection
% Data Using MATLAB by Wail A. Mousa & Abdullatif A. Al-Shuhail
% September 2011.

figure, plot(gx,sx,'s','LineWidth',1)
daspect([1 1 1])
hold on
for i=1:num_shots
    plot(sx(i*num_trace_per_sg(i)),sx(i*num_trace_per_sg(i)),'r*','LineWidth',2)
end
hold off
xlabel('x-axis locations (ft)','FontSize',14)
ylabel('Line number','FontSize',14)
fac=1.5E3;
axis([min(min(sx),min(gx))-fac,max(max(sx),max(gx))+fac,min(min(sx),min(gx))-fac,max(max(sx),max(gx))+fac])
legend('Sources','Receivers',1)
grid
