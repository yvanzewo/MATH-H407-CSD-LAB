close all; clc; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

openinout; %Open the ports of the analog computer.

fc = 40; % cutt off of the butterworth filter
fs = 2*fc; % sampling frequency
Ts = 1/fs; 
Tcycle = Ts;
lengthExp=15; %Set the length of the experiment (in seconds). 
N0=lengthExp/Ts; %Compute the number of points to save the datas.
Data=zeros(N0,2); %Vector saving the datas of velocity and position. If there are several datas to save, change "1" to the number of outputs.
%Data(:,1)store the velocity of motor in V; Data(:,2) store the position in V

% reference position 
%???????????
init_position = 10.8;%cm %choice of the operating point for position
median_position = 13; %cm 
final_position = 15; %cm
%operating point of velocity of motor;
speed_init = 5; %rad/s it is our U_velocity(0) or u(k-1) = u(0) in discrete approximation
%operating point of the current 
i0 = 0.9; % V
%?????????????
prev_speed_command =0;% speed_init;% Should I put it to zero? because is our small signal velocity command which is u(0);

position_ref = init_position*ones(N0,1); % cm ;position reference vector
Data_error = zeros(N0,1);
Current_Commands = zeros(N0,1); %current Vector storing the input sent to the plant(the output of motor controller).
real_position_cm = zeros(N0,1); % position vector storing experiment position in cm not in V
speed_reference = zeros(N0,1); % store the velocity reference compute by the output of the position controller
prev_error = 0; % to store the previous error



kp_m = 2.5;
ki =0.496%1.834 %0.496;
kp_r =0.251%1.31 %0.251;
Dm_s = tf([kp_r ki], [1 0]); % tranfer function of ring controller in continous 
Dm_d = c2d(Dm_s,Ts,'Tustin'); % discrete transfer function of ring controller i
position_cm = 0; % to store the conversion of the voltage position in cm


Xvolt = [3.09,3.27,3.36,3.54,3.73,3.97,4.17,4.52,4.87,5.23,5.74,6.45,7.30,7.82 ];
Xcm = [21.7	20.7	19.7	18.7	17.7	16.7	15.7	14.7	13.7	12.7	11.7	10.7	9.7,8];

% Design a ramp
for j = 1:N0
    if (j<= (N0/3))
         position_ref(j) = init_position + 0.84*(j-1)*Ts;
    else
        position_ref(j) = final_position;
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
    ref_speed = (prev_speed_command + Dm_d.num{1}(1)*error + Dm_d.num{1}(2)*prev_error); %- kp_r*prev_error + (kp_r + Ts*ki)*error; % in rad/s small signal
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
    prev_error = error;
    prev_speed_command = ref_speed;
        
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
