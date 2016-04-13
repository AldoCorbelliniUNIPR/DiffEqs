%% Exercise 1b: analyze the performance of different solvers
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

lambda = theta1+theta2;
options = odeset('Stats','on','Jacobian',-lambda,'RelTol',1e-10);

disp('-----------------------')
disp('ode15s stats:')
tic 
[t,etaAnum]=ode15s(diffeqA,tspan,A0,options);
toc
% etaAana = concentration of A for each element of seqt when the
% differential equation is solved analitically
etaAana=theta2/(theta1+theta2)+ (A0-K)*exp(-(theta1+theta2)*t);
disp('Sum of squares (analytical-numerical)')
disp(sum((etaAnum-etaAana).^2))
disp('Maximum error |analytical-numerical|')
disp(max(abs(etaAnum-etaAana)));


disp('-----------------------')
disp('ode23s stats:')
tic 
[t,etaAnum]=ode23s(diffeqA,tspan,A0,options);
toc
% etaAana = concentration of A for each element of seqt when the
% differential equation is solved analitically
etaAana=theta2/(theta1+theta2)+ (A0-K)*exp(-(theta1+theta2)*t);
disp('Sum of squares (analytical-numerical)')
disp(sum((etaAnum-etaAana).^2))
disp('Maximum error |analytical-numerical|')
disp(max(abs(etaAnum-etaAana)));


disp('-----------------------')
disp('ode23t stats:')
tic 
[t,etaAnum]=ode23t(diffeqA,tspan,A0,options);
toc
% etaAana = concentration of A for each element of seqt when the
% differential equation is solved analitically
etaAana=theta2/(theta1+theta2)+ (A0-K)*exp(-(theta1+theta2)*t);
disp('Sum of squares (analytical-numerical)')
disp(sum((etaAnum-etaAana).^2))
disp('Maximum error |analytical-numerical|')
disp(max(abs(etaAnum-etaAana)));

disp('-----------------------')
disp('ode23tb stats:')
tic 
[t,etaAnum]=ode23tb(diffeqA,tspan,A0,options);
toc
% etaAana = concentration of A for each element of seqt when the
% differential equation is solved analitically
etaAana=theta2/(theta1+theta2)+ (A0-K)*exp(-(theta1+theta2)*t);
disp('Sum of squares (analytical-numerical)')
disp(sum((etaAnum-etaAana).^2))
disp('Maximum error |analytical-numerical|')
disp(max(abs(etaAnum-etaAana)));


disp('-----------------------')
disp('ode45 stats:')
tic 
[t,etaAnum]=ode45(diffeqA,tspan,A0,options);
toc
% etaAana = concentration of A for each element of seqt when the
% differential equation is solved analitically
etaAana=theta2/(theta1+theta2)+ (A0-K)*exp(-(theta1+theta2)*t);
disp('Sum of squares (analytical-numerical)')
disp(sum((etaAnum-etaAana).^2))
disp('Maximum error |analytical-numerical|')
disp(max(abs(etaAnum-etaAana)));


disp('-----------------------')
disp('ode113 stats:')
tic 
[t,etaAnum]=ode113(diffeqA,tspan,A0,options);
toc
% etaAana = concentration of A for each element of seqt when the
% differential equation is solved analitically
etaAana=theta2/(theta1+theta2)+ (A0-K)*exp(-(theta1+theta2)*t);
disp('Sum of squares (analytical-numerical)')
disp(sum((etaAnum-etaAana).^2))
disp('Maximum error |analytical-numerical|')
disp(max(abs(etaAnum-etaAana)));