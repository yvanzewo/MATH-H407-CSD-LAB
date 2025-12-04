close all
fc = 40; % cutt off of the butterworth filter
fs = 2*fc; % sampling frequency
Ts = 1/fs; 
Tcycle = Ts;
lengthExp=20; %Set the length of the experiment (in seconds). 
N0=lengthExp/Ts; %Compute the number of points to save the datas.

init_position = 8;%cm %choice of the operating point for position
median_position = 13; %cm 
final_position = 10.8; %cm
%operating point of velocity of motor;
speed_init = 5; %rad/s it is our U_velocity(0) or u(k-1) = u(0) in discrete approximation
%operating point of the current 
i0 = 0.9; % V

position_ref = init_position*ones(N0,1); % cm ;position reference vector
omega = 0.45*2*pi;

%for ramp and step reference
for j = 1:N0
    if (j<= (N0/4))
         position_ref(j) = init_position + 0.56*(j-1)*Ts;
    elseif (j>N0/4) && (j <= N0/2)
        position_ref(j) = final_position;
    elseif (j>N0/2) && (j<= 3*N0/4)
        position_ref(j) = final_position + 0.44*((j-N0/2)-1)*Ts;
    else
        position_ref(j) = median_position;
    end
end

%for sinus reference control
%  for j = 1:N0
%     if (j<= (N0/4))
%          position_ref(j) = init_position + 1.24*(j-1)*Ts;
%     elseif (j>N0/4) 
%         position_ref(j) = final_position + sin(omega*((j-N0/4)-1)*Ts);
%     end
%  end

time = 0:Tcycle:(N0-1)*Tcycle;

figure 
plot(time,position_ref(:))