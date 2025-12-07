close all; clc;
% Compare experimental results and simulated results for ring control
load('final_ramp_step_trapezoidal_tracking_workspace.mat')

subplot(2,1,1)
plot(time,position_ref(:),'black',time,real_position_cm(:),out.simout.Time,out.simout.Data(:,3),'g')
title('Simulation and experimental results for ring control')
legend("Reference position","Experimental position","Simulated position")
xlabel("Time(s)")
ylabel("Ring position (cm)")
ylim([7.5 20])

subplot(2,1,2)
plot(time,Data(:,1),out.simout.Time,out.simout.Data(:,2))
legend("Experimental Motor velocity","Simulated Motor velocity")
xlabel("Time(s)")
ylabel("Motor velocity (rad/s)")
ylim([2,10])
