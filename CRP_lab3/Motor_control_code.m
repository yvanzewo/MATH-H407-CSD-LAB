close all; clc ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Control System Design Lab: Sample Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

openinout; %Open the ports of the analog computer.

fc = 40; % cutt off of the butterworth filter
fs = 2*fc; % sampling frequency
Ts = 1/fs; 
Tcycle = Ts;
lengthExp=25; %Set the length of the experiment (in seconds). 
N0=lengthExp/Ts; %Compute the number of points to save the datas.
Data=zeros(N0,1); %Vector saving the datas. If there are several datas to save, change "1" to the number of outputs.
speed_ref = 3.5;%V
%choice of the operating point
u0 = 0.9; %y0 will be Data(1,1)
%y0 = 2.1024;% final step
y0 = 2.42; % velocity operating point
Data_ref = y0*ones(N0,1); % speed reference

k_p = 2.5;
DataCommands=ones(N0,1); %Vector storing the input sent to the plant. 
cond=1; %Set the condition variable to 1.
i=1; %Set the counter to 1.
tic %Begins the first strike of the clock.


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while cond==1
    tic
    [in1,in2,in3,in4,in5,in6,in7,in8]=anain; %Acquisition of the measurements.
    % put control law
    Data(i,1) = in1;% store the real velocity output
    if (i*Ts > 5) 
        Data_ref(i) = speed_ref;
    end
    %compute error
     error = k_p*(Data_ref(i) - Data(i,1));
    % send to the plant
    input = error + u0;
    anaout(input,0);
    % store data command(current here) in big signal
    DataCommands(i,1) = input;
        
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
plot(time,Data(:,1),time,DataCommands(:),time,Data_ref(:)); %Plot the experiment (input and output).
legend("Motor velocity","Step current","reference velocity");
