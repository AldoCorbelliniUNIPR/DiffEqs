function NLSminODE(Y,t,theta0)


end

% This file contains the two routines which respectively compute the
% function which  has to be minimized
function dZ=likfx(theta)

Yhat=findyhat(theta,t);
R=Y-Yhat;
dZ=det(R'*R);
end


function Yhat=findyhat(theta,t)
diffeq = @(t,y) [-theta(1)*(y(1)^lambda(1)); theta(1)*(y(1)^lambda(1))-theta(2)*(y(2)^lambda(2)); theta(2)*(y(2)^lambda(2))];

[t,Yhat] = ode45(diffeq,tspan,AB0,options);

end
