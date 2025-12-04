close all; clc;
%PID control with sinus reference

close all; clc; clear;
%Implement PID controller for the ring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

openinout; %Open the ports of the analog computer.

fc = 40; % cutt off of the butterworth filter
fs = 2*fc; % sampling frequency
Ts = 1/fs; 
Tcycle = Ts;
lengthExp=20; %Set the length of the experiment (in seconds). 
N0=lengthExp/Ts; %Compute the number of points to save the datas.
Data=zeros(N0,2); %Vector saving the datas of velocity and position. If there are several datas to save, change "1" to the number of outputs.
%Data(:,1)store the velocity of motor in V; Data(:,2) store the position in V

% initial reference position 
%???????????
init_position = 10.8;%cm %choice of the operating point for position
median_position = 13; %cm 
final_position = 17; %cm
%operating point of velocity of motor;
speed_init = 5; %rad/s it is our U_velocity(0) or u(k-1) = u(0) in discrete approximation
%operating point of the current 
i0 = 0.9; % V

position_ref = init_position*ones(N0,1); % cm ;position reference vector
Data_error = zeros(N0,1);
Current_Commands = zeros(N0,1); %current Vector storing the input sent to the plant(the output of motor controller).
real_position_cm = zeros(N0,1); % position vector storing experiment position in cm not in V
speed_reference = zeros(N0,1); % store the velocity reference compute by the output of the position controller
prev_error1 = 0; % to store error e[k-1]
prev_error2 = 0; %store e[k-2]
prev_speed_command1 =0;%store small velocity ref: u[k-1]
prev_speed_command2 =0;%store small velocity ref: u[k-2]


%Derive the continous and discrete PID controller of Ring
% form is: Dp(s) = kp_r*N_filter*[(ki+s)/s]*[(kd+s)/(N_filter*kd +s)]

kp_m = 4; % 2.5 % controller of motor 
ki = 1.5;
kd = 15.3226; % PID zero  use to cancel a pole of the open loop system
kp_r = 1.77 ; % Gain of PID found using rlocus of open loop system
N_filter = 8; %For the filtered pole Tf = Td/N_filter where Td = 1/kd
part1 = N_filter*kp_r*tf([1 ki],[1 0]); % first part of PID
part2 = tf([1 kd],[1 N_filter*kd]);
Dp_s = part1*part2; % PID in continous
Dp_d = c2d(Dp_s,Ts,'Tustin');% PID in discrete form; Dp(z) = (b0*z²+b1*z+b2)/(z²+ a1*z+a2)
b = Dp_d.num{1}; %[b0 b1 b2]
a = Dp_d.den{1}; % [a0 a1 a2]
position_cm = 0; % to store the conversion of the voltage position in cm


Xvolt = [3.09,3.27,3.36,3.54,3.73,3.97,4.17,4.52,4.87,5.23,5.74,6.45,7.30,7.82 ];
Xcm = [21.7	20.7	19.7	18.7	17.7	16.7	15.7	14.7	13.7	12.7	11.7	10.7	9.7,8];

% Design sinus reference position 
for j = 1:N0
    if (j<= (N0/4))
         position_ref(j) = init_position + 1.24*(j-1)*Ts;
    elseif (j>N0/4) 
        position_ref(j) = final_position + sin(omega*((j-N0/4)-1)*Ts);
    end
 end

cond=1; %Set the condition variable to 1.
i=1; %Set the counter to 1.
%tic %Begins the first strike of the clock.


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while cond==1
    tic
    [in1,in2,in3,in4,in5,in6,in7,in8]=anain; %Acquisition of the measurements.
    % put control law
    Data(i,1) = in1;% store the real velocity output
    Data(i,2) = in2; % store the real position
    speed_rad =in1*2.07; % speed in rad/s not in Volt
    
%     if (i*Ts > 5) % give a step but we will also test with a trapezoidal
%         position_ref(i) = final_position;
%     end
    % convert the voltage position in cm with interpolation
    if (in2 > 7.82)
        position_cm = 8; % in cm
    elseif (in2 < 3.09)
        position_cm = 21.7;
    else
        position_cm = interp1(Xvolt,Xcm,in2,'linear');
    end

    %compute error
    error = position_ref(i,1) - position_cm; % error on the position
    %compute reference velocity to track the motor controller
    % u[k] = -a1*u[k-1] - a2*u[k-2] + b0*e[k] + b1*e[k-1] + b2*e[k-2]
    ref_speed = -a(2)*prev_speed_command1 - a(3)*prev_speed_command2 +b(1)*error ...
                + b(2)*prev_error1 + b(3)*prev_error2;% reference velocity in rad/s small signal; and b(1) is b0 a(2) is a1
    speed_reference(i,1) = ref_speed + speed_init; % (rad/s)reference speed compute by the output of the position controller but has to be in big signal
    if speed_reference(i,1) < 0
        speed_reference(i,1) = 0;
    elseif speed_reference(i,1) > 14
        speed_reference(i,1) = 14; %rad/s
    end
    
    %compute the error of the motor controller: small signal current: i(t)
     error_motor = kp_m*(speed_reference(i,1) - speed_rad)/2.07; % in Volt
    % send to the plant (always big signal)
    input_motor = error_motor + i0;
    if input_motor > 10
        input_motor = 10;
    elseif input_motor < -10
        input_motor = -10;
    end
    anaout(input_motor,0);

    % store data  in big signal
    Current_Commands(i,1) = input_motor; %current to send to the plant
    Data_error(i,1) = error; % the error on the position
    Data(i,1) = speed_rad; % store the real speed of the motor in rad/s
    real_position_cm(i,1) = position_cm; % the real position of ring in cm; Data(:,2) has it but in volt
    %position_ref vector give the tracking ring position
    %???????
    %if initial prev_speed_command is speed_inti so change its initialisation and write
    %speed_reference(i,1) = ref_speed; %in rad/s 

    % store error and speed command for iteration
    prev_error1 = error;
    prev_error2 = prev_error1;
    prev_speed_command1 = ref_speed;
    prev_speed_command2 = prev_speed_command1;
        
    i=i+1;
    t=toc; %Second strike of the clock.
    if t>Tcycle
        disp('Sampling time too small');%Test if the sampling time is too small.
    else
        while toc<=Tcycle %Does nothing until the second strike of the clock reaches the sampling time set.
        end
    end
    if i==N0+1 %Stop condition.
        cond=0;
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

closeinout %Close the ports.
i = i-1;
time = 0:Tcycle:(i-1)*Tcycle;

figure %Open a new window for plot.
plot(time,real_position_cm(:),time,position_ref(:),time,Data(:,1),time,speed_reference); %Plot the experiment (input and output).
legend("Ring position(cm)","reference position(cm)","Motor velocity (rad/s)","Reference velocity(rad/s)");
figure 'Current in Volt sent'
plot(time,Current_Commands(:))
xlabel("time (s)");
ylabel("Current(V)");
