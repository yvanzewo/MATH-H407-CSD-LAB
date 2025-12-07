close all; clc;
%Goal is to design a PID controller for the ring using tool like root locus
%margin.

TF_ring = tf([95.37],[1 18.18 -22.08]);% cm

TF_motor = tf (3.401, [0.9531 1]);

Ts = 1/(2*50);
kp_m = 4;
%kp_m2 = 0.15;
% expression of the inner loop: closed loop of the motor
T_motor = feedback(TF_motor*kp_m,1);
Test_TF = T_motor*TF_ring;

%rlocus(Test_TF);

kp = 1.31;%1;
ki =1.834; %2.21;
k=1 ;%1.404;
Ci = tf([kp ki],[1 0]); % Phase lag controller
%Cd = tf([1 kd],[1 N*kd]); % derivative part of controller
%Dp = N*Cd*Ci; % PID controller without the small k. small k will find with rlocus of open loop 
% form Dp(s) = k* (s + ki)/(s+0.5); where ki> 0.5 
Dp = Ci; % lag controller without the small k. small k will find with rlocus of open loop 
Oloop = k*Dp*Test_TF;
%Oloop.InputDelay = Ts;
Final_closed_loop = feedback(Oloop,1);

figure
%margin(Oloop)
bode(Final_closed_loop)
%step(Final_closed_loop)
%rlocus(Oloop)
