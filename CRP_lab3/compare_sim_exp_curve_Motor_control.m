close all; clc;

load('controlMotor_workspace.mat'); %load data experiment
% data traitment to have velocity output of 2.42V at current input of 0.6V
L = length(out.simout.Data(:,2));
sim_current = zeros(L,1);
exp_current = zeros(N0,1);
sim_current(:,1) =out.simout.Data(:,2);% zeros(N0,1);
exp_current(:,1) = DataCommands(1:N0,1);%zeros(N0,1);

N = (5/Ts);
N3= (5.012/0.001) +200;
N2 =N;

sim_current(1:N3,1) = out.simout.Data(1:N3,2) - 0.3;

exp_current(1:N2,1) = DataCommands(1:N2,1) - 0.3;

% overlay the simulation curve and the real experiment plant for the motor
subplot(2,1,1) 
plot(out.simout.Time,out.simout.Data(:,1),time,Data(:,1),time,Data_ref(:))
%plot(out.simout.Time,out.simout.Data(:,1),time,Data(:,1),out.simout.Time,out.simout.Data(:,2),time,DataCommands(:),time,Data_ref(:))
title('Validation of the motor control with a P controller of kp=2.5');
legend("Simulated velocity"," experimental velocity","Reference velocity")
ylim([1.5,3.8]);
xlabel("Time(s)");
ylabel("Motor velocity(V)");

subplot(2,1,2) % for current experimental and simulated
plot(out.simout.Time,sim_current(:),time,exp_current(:))
legend("Simulated current"," experimental Current")
ylim([0,4])
xlabel("Time(s)");
ylabel("Motor current(V)");