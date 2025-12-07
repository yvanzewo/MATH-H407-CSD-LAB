close all; clc;
fc = 40; % cutt off of the butterworth filter
fs = 2*(fc+10); % sampling frequency
Ts = 1/fs; 
Tcycle = Ts;
lengthExp=55; %Set the length of the experiment (in seconds). 
N0=lengthExp/Ts; %Compute the number of points to save the datas.

init_position = 8;%cm %choice of the operating point for position
final_position = 17; %cm 
final_position2 = 13; %cm 
eq_position = 10.8; %cm equilibrium point
%operating point of velocity of motor;
speed_init = 5; %rad/s it is our U_velocity(0) or u(k-1) = u(0) in discrete approximation
%operating point of the current 
i0 = 0.9; % V

position_ref = init_position*ones(N0,1); % cm ;position reference vector
omega = 0.45*2*pi;
% time to design the reference
% time to design the reference
T1 = 5; %second
N1 = (T1/lengthExp)*N0;
T2 = 20;
N2 = (T2/lengthExp)*N0;
T3 = 25;
N3 = (T3/lengthExp)*N0;
T4 = 40;
N4 = (T4/lengthExp)*N0;
T5 = 45;
N5 = (T5/lengthExp)*N0;
%ref to reach equilibrium point
% for j = 1:N0
%     if (j<= N1)
%          position_ref(j) = init_position + ((eq_position-init_position)/(T1))*(j-1)*Ts;
%     else
%         position_ref(j) = eq_position;
%     end
% end 

    
% ref with multiple ramp and step
% for j = 1:N0
%     if (j<= N1)
%          position_ref(j) = init_position + ((eq_position-init_position)/(T1))*(j-1)*Ts;
%     elseif (j > N1) && (j <= N2)
%         position_ref(j) = eq_position;
%     elseif (j>N2) && (j<= N3)
%         position_ref(j) = eq_position + ((final_position-eq_position)/(T3-T2))*((j-N2)-1)*Ts;
%     elseif (j>N3) && (j<= N4) 
%         position_ref(j) = final_position;
%     elseif (j>N4) && (j<= N5) 
%         position_ref(j) = final_position + ((final_position2-final_position)/(T5-T4))*((j-N4)-1)*Ts;
%     else
%         position_ref(j) = final_position2;
%     end
% end
% sinus reference design
for j = 1:N0
    if (j<= N1)
         position_ref(j) = init_position + ((eq_position-init_position)/(T1))*(j-1)*Ts;
    elseif (j > N1) && (j <= N2)
        position_ref(j) = eq_position;
    elseif (j>N2) && (j<= N3)
        position_ref(j) = eq_position + ((final_position-eq_position)/(T3-T2))*((j-N2)-1)*Ts;
    elseif (j>N3) && (j<= N4) 
        position_ref(j) = final_position;
    else
        
         position_ref(j) = final_position + sin(omega*((j-N4)-1)*Ts);
    end
end

%for ramp and step reference
% for j = 1:N0
%     if (j<= (N0/4))
%          position_ref(j) = init_position + 0.56*(j-1)*Ts;
%     elseif (j>N0/4) && (j <= N0/2)
%         position_ref(j) = final_position;
%     elseif (j>N0/2) && (j<= 3*N0/4)
%         position_ref(j) = final_position + 0.44*((j-N0/2)-1)*Ts;
%     else
%         position_ref(j) = median_position;
%     end
% end

%for sinus reference control
%  for j = 1:N0
%     if (j<= (N0/4))
%          position_ref(j) = init_position + 1.24*(j-1)*Ts;
%     elseif (j>N0/4) 
%         position_ref(j) = eq_position + sin(omega*((j-N0/4)-1)*Ts);
%     end
%  end

time = 0:Tcycle:(N0-1)*Tcycle;

figure 
plot(time,position_ref(:))