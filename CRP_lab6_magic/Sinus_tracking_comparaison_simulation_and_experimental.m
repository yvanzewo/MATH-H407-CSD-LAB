close all; clc;
% compare sinus reference simulation and experimental 
load('sinus_freq_0.07_amp_1_workspace.mat')

Ts_sim = 0.01;
exp_start = 40/Ts;
exp_end = 75/Ts;
sim_start = 42.9/Ts_sim;
sim_end = (42.9 +35)/Ts_sim;


interesting_exp_time = time(exp_start:exp_end);
interesting_exp_time = interesting_exp_time -40;
interesting_sim_time = out.simout.Time(sim_start:sim_end);
interesting_sim_time = interesting_sim_time -42.9;

figure
plot(interesting_exp_time,real_position_cm(exp_start:exp_end,1),interesting_sim_time,out.simout.Data(sim_start:sim_end,3),interesting_exp_time,position_ref(exp_start:exp_end,1))
legend("Experimental position","Simulated position","Reference position");
title('Sinus reference tracking f = 0.07Hz')
xlabel ('Time(s)')
ylabel('Position (cm)')
xlim ([0 35])
