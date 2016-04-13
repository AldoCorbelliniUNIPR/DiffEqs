%% exercise 1c: incorporate the mass balanc contrainst [A]+[B]=M

close all
% time range
tspan = [0 5];

% In this exercise theta1 and theta2 are known constants
theta1=1;
theta2=2;
K=theta2/(theta1+theta2);

% Differential equation for A
% d[A]/dt = -theta1 [A] + theta2 [B]
% Differential equation for B
% d[B]/dt =  theta1 [A] - theta2 [B]
% Constraint 0=[A]+[B]-1

% Remark: an anonymous function is created
% The output of the anonymous function is of size 3-by-1
usingMassMatrix=1;
if usingMassMatrix==1
diffeqAB = @(t,y) [-theta1*y(1)+theta2*y(2); theta1*y(1)-theta2*y(2); y(1)+y(2)-1];
options = odeset('Mass',M,'RelTol',1e-4,'AbsTol',[1e-6 1e-6 1e-6]);
M = [1 0 0; 0 1 0; 0 0 0];
AB0=[1; 0; 0];

else
 diffeqAB = @(t,y) [-theta1*y(1)+theta2*y(2); theta1*y(1)-theta2*y(2)];
 options = odeset('RelTol',1e-4,'AbsTol',[1e-6 1e-6]);
 AB0=[0.8; 0.2];
end

% A0 = concentration of A at t=0
% AB0=[2.1; 1];
%AB0=[3; 2];
% Set the tolerance.
% Remark: please note that the length of vector increase as RElTol
% decreases
% options = odeset('RelTol',1e-13,'stats','on');

% etaABnum= in this case etaABnum is a matrix with three columns
% 1st column = concentration of A at the integration points t 
% 2nd column = concentration of B at the integration points t 
% 3rd colum  = a vector of zeros (reflecting the constraint)
[t,etaABnum] = ode45(diffeqAB,tspan,AB0,options);

