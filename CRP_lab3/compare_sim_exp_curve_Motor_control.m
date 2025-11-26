close all; clc;



% overlay the simulation curve and the real experiment plant for the motor
figure

plot(out.simout.Time,out.simout.Data(:,1),time,Data(:,1),out.simout.Time,out.simout.Data(:,2),time,DataCommands(:),time,Data_ref(:))
title('Validation of the motor control with a P controller of kp=2.5');
legend("Simulated velocity","velocity experiment","Simulated current","Current experiment","Reference velocity")
xlabel("Time(s)");
ylabel("Amplitude(V)");
