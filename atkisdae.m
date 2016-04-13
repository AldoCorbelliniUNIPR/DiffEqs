function res = atkisdae(t,y)
theta=[1 2];
res = [-theta(1)*y(1) + theta(2)*y(2);
   theta(1)*y(1) - theta(2)*y(2);
   y(1) + y(2) - 1];
end
