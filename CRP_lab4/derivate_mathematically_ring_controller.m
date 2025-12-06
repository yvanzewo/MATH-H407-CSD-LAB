close all; clc;
%Goal is to design a PID controller for the ring using tool like root locus
%margin.

TF_ring = tf([95.37],[1 18.18 -22.08]);% cm
%TF_ring2 = tf(0.97,[1 12.55 -41.94]);
TF_motor = tf (3.401, [0.9531 1]);
%TF_motor2 = tf(729.44,[1.31 1]);
Ts = 1/(2*50);
kp_m = 4;
%kp_m2 = 0.15;
% expression of the inner loop: closed loop of the motor
T_motor = feedback(TF_motor*kp_m,1);
Test_TF = T_motor*TF_ring;

%rlocus(Test_TF);
kd = T_motor.den{1}(2)/T_motor.den{1}(1)
kp = 1;%1;
ki =1.4; %2.21;
k=1 ;%1.404;
%N = 8;
Ci = tf([kp ki],[1 0.5]); % Phase lag controller
Cd = tf([1 kd],[1 N*kd]); % derivative part of controller
%Dp = N*Cd*Ci; % PID controller without the small k. small k will find with rlocus of open loop 
% part1 = N*850.802*Ci;% first part of openloop ring controller
% part2 = tf([1],[1 N*kd]); % second
% part3 = tf([1],[1 kd]); %third
% part4 = tf([1],[1 -1.14]);
%Oloop = Dp*Test_TF; %open loop of ring controller
%Oloop = part1*part2*part3*part4;
%Oloop.InputDelay = Ts;
Dp = Ci; % PID controller without the small k. small k will find with rlocus of open loop 
Oloop = k*Dp*Test_TF;
%Oloop.InputDelay = Ts;
Final_closed_loop = feedback(Oloop,1);

figure
%margin(Oloop)
%bode(Final_closed_loop)
%step(Final_closed_loop)
rlocus(Oloop)
