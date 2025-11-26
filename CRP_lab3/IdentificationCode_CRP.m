u=DataCommands.';
y=Data(:,1).';
u_offset=u0.'; % input Operating point
y_offset = Data(1,1).'; % output operating point
SystemOrder=[0 1]; %Number of zeros and of poles (0 and 1), respectively.
sysIdent=IdentifySystem_CRP(u-u_offset,y-y_offset,SystemOrder,Ts)
figure
plot(time,y-y_offset,'.');
hold on;
lsim(sysIdent,u-u_offset,time);


% y_offset = Data(20/Ts.',1).'; % output operating point
% SystemOrder=[0 1]; %Number of zeros and of poles (0 and 1), respectively.
% sysIdent=IdentifySystem_CRP(u-u_offset,y-y_offset,SystemOrder,Ts);
% figure
% plot(time((20/Ts.'):(lengthExp.'/Ts.'),y-y_offset,'.');
% hold on;
% lsim(sysIdent,u-u_offset,time(20/Ts.':lengthExp.'/Ts.');