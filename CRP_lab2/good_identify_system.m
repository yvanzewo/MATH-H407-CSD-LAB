close all; clear; clc;
load('data for tf of motor 2.mat'); %load experiment data
u=DataCommands.';
y=Data(:,1).';
start = 19/Ts; %starting time we want to record our linearized transfert fct
finish = 30/Ts;%Final time to record the Transfer function
u_offset=u(1,start); % input Operating point (normally is u_offset = u0 = 0.48 V choosen point to have 5 rad/s)
y_offset = y(1,start); %y0 output operating point
SystemOrder=[0 1]; %Number of zeros and of poles (0 and 1), respectively.
sysIdent=IdentifySystem_CRP(u(start:finish)-u_offset,y(start:finish)-y_offset,SystemOrder,Ts)
figure
plot(time(start:finish),y(start:finish)-y_offset,'.');
hold on;
lsim(sysIdent,u(start:finish)-u_offset,time(start:finish));

%Note that our experiment in lab said for an input current of 0.48V, will
%have a velocity of around 2.42V :the value obtained after convertion that states 1V
%is 2.07 rad/s so for 5 rad/s is 2.42V so operation point is (u0,y0) =
%(0.48V,2.42V)
% But workspace shown that for y_offset = 2.1V so the operating point is 
%(u0,y0) = (0.48V,2.1V)