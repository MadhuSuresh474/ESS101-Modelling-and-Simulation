clc;                                     
clear all;
%% 
%Inputs for the linspace 
a=[0.1,0.01,0.001,0.0001];
dt=0.1 ;                    %step size  
dt2=0.01;
dt3=0.001;
dt4=0.0001;

% 10^(-1) being the magnitude of the timesteps will take us 21 time steps to reach 2 seconds from  seconds
t = linspace(0,2,2/dt+1); 
t2 = linspace(0,2,2/dt2+1)
t3 = linspace(0,2,2/dt3+1)
t4 = linspace(0,2,2/dt4+1)

% creating an empty matrix to store our resultant RK estimated resultant points.
y1=zeros(2/dt+1,1); 
y2=zeros(2/dt2+1,1);
y3=zeros(2/dt3+1,1);
y4=zeros(2/dt4+1,1);
        
%Function obeying xdot=lamda.x hence input is not considered
lam= -2;                
f_x = @(t) lam*exp(lam*t);
f_x2 = @(t2) lam*exp(lam*t2);
f_x3 = @(t3) lam*exp(lam*t3);
f_x4 = @(t4) lam*exp(lam*t4);

%Actual function vector
f_act= exp(lam*t);
f_act2= exp(lam*t2);
f_act3= exp(lam*t3);
f_act4= exp(lam*t4);
%f_act5= exp(lam*t5);
% Initial point of the function to compute from.
y(1) = 1;          
y2(1) = 1;
y3(1) = 1;
y4(1) = 1;

%%
%Considering RK4 for the ODE 
yRK4 = rangekutta4mod(dt,y,t,f_x);   %Calling the rangekutta 4 function to estimate the function f_x
localerrorvector4=yRK4-f_act;      %Localerror at each step of the itration in vecor form;
globalerror4=sum(localerrorvector4) %Globalerror at each step of the itration in vecor form;
%%
%Considering RK4 for the ODE dt2
yRK42 = rangekutta4mod(dt2,y2,t2,f_x2);   %Calling the rangekutta 4 function to estimate the function f_x
localerrorvector42=yRK42-f_act2.';      %Localerror at each step of the itration in vecor form;
globalerror42=sum(localerrorvector42) %Globalerror at each step of the itration in vecor form;
%%
%Considering RK4 for the ODE dt2
yRK43 = rangekutta4mod(dt3,y3,t3,f_x3);   %Calling the rangekutta 4 function to estimate the function f_x
localerrorvector43=yRK43-f_act3.';      %Localerror at each step of the itration in vecor form;
globalerror43=sum(localerrorvector43) %Globalerror at each step of the itration in vecor form;
%%
%Considering RK4 for the ODE dt2
yRK44 = rangekutta4mod(dt4,y4,t4,f_x4);   %Calling the rangekutta 4 function to estimate the function f_x
localerrorvector44=yRK44-f_act4.';      %Localerror at each step of the itration in vecor form;
globalerror44=sum(localerrorvector44) %Globalerror at each step of the itration in vecor form;

%%
%Considering RK2 for the ODE 
yRK2 = rangekutta2(dt,y,t,f_x);      %Calling the rangekutta 4 function to estimate the function f_x
localerrorvector2=yRK2-f_act;      %Localerror at each step of the itration in vecor form;
globalerror2=sum(localerrorvector2) %Globalerror at each step of the itration in vecor form;
%%
%Considering RK2 for the ODE dt2
yRK22 = rangekutta2(dt2,y2,t2,f_x2);      %Calling the rangekutta 4 function to estimate the function f_x
localerrorvector22=yRK22-f_act2.';      %Localerror at each step of the itration in vecor form;
globalerror22=sum(localerrorvector22) %Globalerror at each step of the itration in vecor form;
%%
%Considering RK2 for the ODE dt2
yRK23 = rangekutta2(dt3,y3,t3,f_x3);      %Calling the rangekutta 4 function to estimate the function f_x
localerrorvector23=yRK23-f_act3.';      %Localerror at each step of the itration in vecor form;
globalerror23=sum(localerrorvector23) %Globalerror at each step of the itration in vecor form;
%%
%Considering RK2 for the ODE dt2
yRK24 = rangekutta2(dt4,y4,t4,f_x4);      %Calling the rangekutta 4 function to estimate the function f_x
localerrorvector24=yRK24-f_act4.';      %Localerror at each step of the itration in vecor form;
globalerror24=sum(localerrorvector24) %Globalerror at each step of the itration in vecor form;

%%
%Considering RK1 or explicit euler scheme for the ODE 
yRK1 = rangekutta1euler(dt,y,t,f_x); %Calling the rangekutta 4 function to estimate the function f_x
localerrorvector1=yRK1-f_act;      %Localerror at each step of the itration in vecor form;
globalerror1=sum(localerrorvector1) %Globalerror at each step of the itration in vecor form;
%%
%Considering RK1 or explicit euler scheme for the ODE dt2
yRK12 = rangekutta1euler(dt2,y2,t2,f_x2); %Calling the rangekutta 4 function to estimate the function f_x
localerrorvector12=yRK12-f_act2.';      %Localerror at each step of the itration in vecor form;
globalerror12=sum(localerrorvector12) %Globalerror at each step of the itration in vecor form;
%%
%Considering RK1 or explicit euler scheme for the ODE dt2
yRK13 = rangekutta1euler(dt3,y3,t3,f_x3); %Calling the rangekutta 4 function to estimate the function f_x
localerrorvector13=yRK13-f_act3.';      %Localerror at each step of the itration in vecor form;
globalerror13=sum(localerrorvector13) %Globalerror at each step of the itration in vecor form;
%%
%Considering RK1 or explicit euler scheme for the ODE dt2
yRK14 = rangekutta1euler(dt4,y4,t4,f_x4); %Calling the rangekutta 4 function to estimate the function f_x
localerrorvector14=yRK14-f_act4.';      %Localerror at each step of the itration in vecor form;
globalerror14=sum(localerrorvector14) %Globalerror at each step of the itration in vecor form;

%%
%plots fot dt and its RK order
figure(1);
plot(t,f_act,'-o',t,yRK4,'-o',t,yRK2,'-o',t,yRK1,'-o')
title('Actual plot vs RK estimates');
xlabel('Time t');
ylabel('Solution y');
legend('fact','yRK4','yRK2','yRK1')
GBE4=[abs(globalerror4),abs(globalerror42),abs(globalerror43),abs(globalerror44)]  
GBE2=[abs(globalerror2),abs(globalerror22),abs(globalerror23),abs(globalerror24)]
GBE1=[abs(globalerror1),abs(globalerror12),abs(globalerror13),abs(globalerror14)]
figure(2);
loglog(a,GBE4,'-o',a,GBE2,'-o',a,GBE1,'-o')
title('Error vs deltaT');
xlabel('Time t');
ylabel('Solution y');
legend('GBE4','GBE2','GBE1')
% Thus from figure (2) we can conclude that accuracy increases drastically
% with small variation in delta T in RK4 (higher order) compared to RK(1)
% where the accuracy chenges mildly even with much lower valves of delta T.