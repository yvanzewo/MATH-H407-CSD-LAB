close all; clear; clc;
% Acess the stability margin of the control motor simulation;

load('data for tf of motor 2.mat'); %load exp data
k_p = 2.5; % proportional gain for P controller
D_motor = k_p; 
TF_motor = tf([3.401], [0.9531 1]); % Transfer function of Motor around the choosen operating point
Delay = Ts; %  Delay define the calculation time of ADC and DAC; also the ZOH impact
%Anti alaising filter: 2nd order Butterworth filter
omega_n = 2*pi*fc ;
zeta = 0.7;
Filter_AA = tf([omega_n^2],[1 2*omega_n*zeta omega_n^2]); % TF of filter
OP_motor_1 = D_motor*TF_motor; % open loop of motor without delay and filter;
OP_motor_1.InputDelay = Delay; % Add delay on the inner loop of motor
OP_motor_2 = OP_motor_1*Filter_AA; % open loop of motor with filter and delay

[Gm1,Pm1] = margin(OP_motor_1); % gain and phase margin without filter
[Gm2,Pm2] = margin (OP_motor_2); % gain and phase margin with filter and delay
figure;
margin(OP_motor_1);
figure;
margin (OP_motor_2);