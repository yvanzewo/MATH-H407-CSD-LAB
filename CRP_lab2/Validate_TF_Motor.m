close all; clc;% don't put clear because it clear the workspace

% Validation of transfer function of motor by overlay the simulated output
% with real output but in small signal because the TF obtained it is in
% small signal;

load('data for tf of motor 3.mat'); %load experiment data


% y=Data(:,1).';
% start = 19/Ts; %starting time we want to record our linearized transfert fct
% finish = 30/Ts;%Final time to record the Transfer function
% y_offset = y(1,start); %y0 output operating point

y_small_signal = y(start:finish)-y_offset; % output small signal that interested us 
time_interested = time(start:finish)-start*Ts;

figure
plot(out.simout.Time ,out.simout.Data(:,1),time_interested,y_small_signal)
title('Transfer function validation of motor');
legend("Simulated function","experiment function")
xlabel("Time(s)");
ylabel("Amplitude(V)");
xlim([0 11]);


