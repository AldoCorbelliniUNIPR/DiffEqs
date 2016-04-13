%% exercise 1a: compare numerical and analytical solution 
close all
% time range
tspan = [0 5];

% In this exercise theta1 and theta2 are known constants
theta1=1;
theta2=2;
K=theta2/(theta1+theta2);

% Differential equation for A
% d[A]/dt = -theta1 [A] + theta2 (1-A) = -(theta1+theta2)*A
% Remark: an anonymous function is created
diffeqA = @(t,A) -theta1*A+theta2*(1-A);

% A0 = concentration of A at t=0
A0=1;
% Set the tolerance.
% Remark: please note that the length of vector increase as RElTol
% decreases
options = odeset('RelTol',1e-3,'stats','on');

% etaAnum= concentration of A at the integration points t using ode45
[t,etaAnum] = ode45(diffeqA,tspan,A0,options);

seqt=t;

% etaAana = concentration of A for each element of seqt when the
% differential equation is solved analitically
etaAana=theta2/(theta1+theta2)+ (A0-K)*exp(-(theta1+theta2)*seqt);


% Compare the two solutions
subplot(2,1,1)
% Results from the analytical solution
plot(seqt,etaAana,'-x')
title('Differential equation for A: analytical solution')

subplot(2,1,2)
% Results from the numerical solution
tspan = [0 5];
plot(t,etaAnum,'-x')
xlabel('t')
title('Differential equation for A: numerical solution')
if prin ==1
 print -depsc figs/ex1a.eps;
end
%% Differential equation for B
close all
% d[B]/dt = -(theta1+theta2)*B +theta1
% Remark: an anonymous function is created
diffeqB = @(t,B) -(theta1+theta2)*B+theta1;

% B0 = concentration of B at t=0
B0=1;
% seqt = sequence of time points
seqt=0:0.05:5;

% etaBana = concentration of A for each element of seqt when the
% differential equation is solved analitically
etaBana=1-theta2/(theta1+theta2)+ (B0-K)*exp(-(theta1+theta2)*seqt);

% etaAnum= concentration of A at the integration points t using ode45
[t,etaBnum] = ode45(diffeqB,tspan,B0);

% Compare the two solutions
subplot(2,1,1)
% Results from the analytical solution
plot(seqt,etaBana,'-x')
title('Analytical solution')

subplot(2,1,2)
% Results from the numerical solution
tspan = [0 5];
plot(t,etaBnum,'-x')
xlabel('t')
title('Numerical solutions')


