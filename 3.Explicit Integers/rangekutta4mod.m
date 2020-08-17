function [y]= rangekutta4mod(dt,y,t,f_x)
k = zeros(1,4);
      
A1=[0 0 0 0;0.5 0 0 0;0,0.5,0,0;0 0 1 0];
B1=[1/6 1/3 1/3 1/6];
C1=[0 0.5 0.5 1];
        
for i=1:(length(t)-1)  
    
    k(1) = f_x(t(i));
    k(2) = f_x(t(i)+A1(2,1)*dt);
    k(3) = f_x((t(i)+A1(3,2)*dt));
    k(4) = f_x((t(i)+A1(4,3)*dt));
    
    y(i+1) = y(i) + B1*k.'*dt;  % main equation
end