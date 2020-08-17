clc;                                     
clear all;
%% 
%Inputs for the linspace 
%a=[0.1,0.01,0.001,0.0001];
dt=0.01 ;                    %step size  

% 10^(-1) being the magnitude of the timesteps will take us 21 time steps to reach 2 seconds from  seconds
t = linspace(0,2,2/dt+1); 

% creating an empty matrix to store our resultant RK estimated resultant points.
%y1=zeros(2/dt+1,1); 

%Function obeying xdot=lamda.x hence input is not considered
%lam= -2;                
%f_x = @(t) (-2)*exp((-2)*t);
f_x = @(x,y)[y ; 5*(1-x^2)*y-x]
%f_act= exp(-2*t);


y(1) = 1;          
x(1) = 1;
%function [y]= rangekutta2mod(dt,y,t,f_x)
syms k1 k2;
k1 = sym('k1');
k2 = sym('k2');
     
A2=[0.25 0.25-sqrt(3)/6;0.25+sqrt(3)/6 0.25];
B2=[0.5-sqrt(3)/6 0.5-sqrt(3)/6];
C2=[0.5 0.5];

 
for i=1:(length(t)-1)
    k1(1)=y(i);
    k2(1)=y(i);
    j=1;
    
    r = sym('r');
    
    
  
    
    m=i;
    u=1;
    while u>1e-2
        
    
        r = [k1(m)-f_x(x(i)+dt*(A2(1,1)*k1(m)+A2(1,2)*k2(m)),y(i)+dt*(A2(1,1)*k1(m)+A2(1,2)*k2(m))),k2(m)-f_x(x(i)+dt*(A2(2,1)*k1(m)+A2(2,2)*k2(m)),y(i)+dt*(A2(2,1)*k1(m)+A2(2,2)*k2(m)))]
        d0=y(i);
        d1=x(i);
        d2=k1(m)
        d3=k2(m)
        k=[k1(m),k2(m)];
        b=[k1;k2]
       
        syms k11 k22
        c=[k11,k22]
        syms F pp qq;
        
        %F= jacobian([k11-f_x(pp+dt*(A2(1,1)*k11+A2(1,2)*k22)),k22-f_x(pp+dt*(A2(2,1)*k11+A2(2,2)*k22))],c);
        %F= jacobian([k11-f_x(pp+dt*(A2(1,1)*k11+A2(1,2)*k22),qq+dt*(A2(1,1)*k11+A2(1,2)*k22)),k22-f_x(pp+dt*(A2(2,1)*k11+A2(2,2)*k22),qq+dt*(A2(2,1)*k11+A2(2,2)*k22))],c)
        %F= jacobian([k11-f_x(pp+dt*(A2(1,1)*k11+A2(1,2)*k22),qq+dt*(A2(1,1)*k11+A2(1,2)*k22)) , k22-f_x(pp+dt*(A2(2,1)*k11+A2(2,2)*k22),qq+dt*(A2(2,1)*k11+A2(2,2)*k22))],c)
        F= jacobian([(k11-(f_x((pp+dt*(A2(1,1)*k11+A2(1,2)*k22)),(qq+dt*(A2(1,1)*k11+A2(1,2)*k22))))) ; (k22-(f_x((pp+dt*(A2(2,1)*k11+A2(2,2)*k22)),(qq+dt*(A2(2,1)*k11+A2(2,2)*k22)))))],c)
        %F=jacobian([k11-f_x(x+dt*(A2(1,1)*k11+A2(1,2)*k22),y+dt*(A2(1,1)*k11+A2(1,2)*k22)),k22-f_x(x+dt*(A2(2,1)*k11+A2(2,2)*k22),y+dt*(A2(2,1)*k11+A2(2,2)*k22))],c)
        %syms a b c d e
        %M = [a b c; d e b; e c a];
        %Ms = subs(M, {a b c d e}, {1 3 5 7 9});
        FE= double(subs(F,{pp,qq,k11,k22},{d1,d0,d2,d3}))
        %y = feval(F,x1,x2)
        dK=sym('dK');
        %dK=double(-(inv(FE)*[k1(m)-f_x(t(i)+dt*(A2(1,1)*k1(m)+A2(1,2)*k2(m)));k2(m)-f_x(t(i)+dt*(A2(2,1)*k1(m)+A2(2,2)*k2(m)))]));
        dK=double(-(inv(FE)*[k1(m)-f_x(x(i)+dt*(A2(1,1)*k1(m)+A2(1,2)*k2(m)),y(i)+dt*(A2(1,1)*k1(m)+A2(1,2)*k2(m)));k2(m)-f_x(x(i)+dt*(A2(2,1)*k1(m)+A2(2,2)*k2(m)),y(i)+dt*(A2(2,1)*k1(m)+A2(2,2)*k2(m)))]))
        %S= solve(eqn);
        j=m+1;
        d4= double(dK.'+ [d2,d3])
        k1(j)= d4(1)
        k2(j)= d4(2)
        m=m+1;
        %r = [k1(m)-f_x(t(i)+dt*(A2(1,1)*k1(m)+A2(1,2)*k2(m)));k2(m)-f_x(t(i)+dt*(A2(2,1)*k1(m)+A2(2,2)*k2(m)))];
        r1 = [k1(m)-f_x(x(i)+dt*(A2(1,1)*k1(m)+A2(1,2)*k2(m)),y(i)+dt*(A2(1,1)*k1(m)+A2(1,2)*k2(m))),k2(m)-f_x(x(i)+dt*(A2(2,1)*k1(m)+A2(2,2)*k2(m)),y(i)+dt*(A2(2,1)*k1(m)+A2(2,2)*k2(m)))]
        u=abs(min(r1));
    end
    y(i+1)=y(i)+dt*0.5*k1(m)+dt*0.5*k2(m);
    x(i+1)=x(i)+dt*0.5*k1(m)+dt*0.5*k2(m);
end
%plots fot dt and its RK order
figure(1);
plot(t,y,'-o',t,x,'-o');
title('Actual plot vs RK estimates');
xlabel('Time t');
ylabel('Solution y');
%legend('fact','yRK4','yRK2','yRK1'