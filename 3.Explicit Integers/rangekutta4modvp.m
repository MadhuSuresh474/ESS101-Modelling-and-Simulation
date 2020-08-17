function [x,y]= rangekutta4modvp(dt,x,y,g,t,df_x_y)
k = zeros(2,4);
      
A1=[0 0 0 0;0.5 0 0 0;0,0.5,0,0;0 0 1 0];
B1=[1/6 1/3 1/3 1/6];
C1=[0 0.5 0.5 1];
%A1=[0 0 0 0;0.9 0 0 0;0,0.9,0,0;0 0 0.8 0];
%B1=[0.9 0.1 0.1 0.9];
%C1=[0 0.2 0.2 3];
for i=1:(length(t)-1)
     
     k(:,1) = df_x_y(x(i),y(i));
     k(:,2) = df_x_y(x(i)+A1(2,1)*dt,y(i)+A1(2,1)*dt);
     k(:,3) = df_x_y((x(i)+A1(3,2)*dt),(y(i)+A1(3,2)*dt));
     k(:,4) = df_x_y((x(i)+A1(4,3)*dt),(y(i)+A1(4,3)*dt));
     x(i+1) = x(i) + B1*k(1,:).'*dt;
     y(i+1) = y(i) + B1*k(2,:).'*dt;
     %g(1,i+1) = g(1,i) + B1*k(1,:).'*dt;
     %g(2,i+1) = g(1,i) + B1*k(1,:).'*dt;% main equation
end