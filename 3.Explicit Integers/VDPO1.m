function f = VDPO1(t,in2)
%VDPO1
%    F = VDPO1(T,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    18-Oct-2019 23:13:28

y1 = in2(1,:);
y2 = in2(2,:);
f = [y2;-y1-y2.*(y1.^2.*5.0-5.0)];