clc;
clear;
time = [0 25];
y0 = [1; 0];
u=5;
ode = @(t,y) vanderpol(t,y,u);
[t,y] = ode45(ode, [0 25], [2; 0]);
figure(1)
plot(t,y);
xlabel('Time t');
ylabel('Solution y');
legend('x','y');
title('Solution of van der Pol Equation from ODE45');
%%
dt=0.01 ;                    %step size  
t = linspace(0,25,25/dt+1); 

x=zeros(1,25/dt+1);   %y(1)=x     
y=zeros(1,25/dt+1);  %y(2)=y   

x(1) = 2;
y(1) = 0;  
g=[x;y];

df_x_y =  zeros(2,25/dt+1); 
df_x_y = @(x,y)[y; 5*(1-x^2)*y-x]


%%
%Considering RK4 for the ODE f_x
[xRK4,yRK4] = rangekutta4modvp(dt,x,y,g,t,df_x_y);   %Calling the rangekutta 4 function to estimate the vanderpol functiobn

figure(2);
plot(t,yRK4,'-o',t,xRK4,'-o')
title('Solution of van der Pol Equation from RK');
xlabel('Time t');
ylabel('Solution y');
legend('y','x')

% consideting dt=0.001 yeilds us an RK4 estimate equal to its ode45
% estimate for the vanderpol function.