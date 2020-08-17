function [y]= rangekutta2(dt,y,t,f_x)
k = zeros(1,2);
A2=[0 0;0.5 0];
B2=[0 1];
C2=[0 0.5];
        
for i=1:(length(t)-1)  
    
    k(1) = f_x(t(i));
    k(2) = f_x(t(i)+A2(2,1)*dt);
          

    y(i+1) = y(i) + B2*k.'*dt; 
end