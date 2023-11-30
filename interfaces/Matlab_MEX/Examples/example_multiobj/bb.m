function [eval ] = bb(x)

obj1=0;
obj2=0;
for i=1:2
  obj1 = obj1 -10*exp(-0.2*(x(i)*x(i)+x(i+1)*x(i+1))^0.5);
end

for i=1:3
   obj2 = obj2 + abs(x(i))^0.8+5.0*sin(x(i)^3.0);
end

eval=[obj1 , obj2];
