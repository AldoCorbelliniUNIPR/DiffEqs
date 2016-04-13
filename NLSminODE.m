function thetahat=NLSminODE(Y,t,theta0,AB0,lambda,options)

loglik=@likfx;
% fminsearch is the minimization routine
% loglik will contain the function which will have to be minimized
% theta0 is the vector of initial parameters
[thetahat,fval,exitflag]  = fminsearch(loglik,theta0);


% This file contains the two routines which respectively compute the
% function which  has to be minimized
    function dZ=likfx(theta)
        % likfx sets up the target function
        
        % Parameter transformation to help convergence
        % TODISCUSS
        % theta=normcdf(theta);
        theta=sin(theta);
        % disp(theta)
        
        % Call routines which find fitted values in correspondence of t
        % for vector thetat and starting conditions AB0
        Yhat=findyhat(theta,t,AB0);
        % R = Matrix of residuals
        R=Y-Yhat;
        % dZ = determinant of R (the value of dZ must be minimized)
        dZ=det(R'*R);
        disp('Det')
        disp(dZ)
    end


    function Yhat=findyhat(theta,t,AB0)
        % Find fitted values solving ODE
        
        % Define the ODE as function of elements of vector theta
        % Note that elements of theta change at each iteration
         diffeq = @(t,y) [-theta(1)*(y(1)^lambda(1)); theta(1)*(y(1)^lambda(1))-theta(2)*(y(2)^lambda(2)); theta(2)*(y(2)^lambda(2))];
        
         % Given the elements of theta solve the ODE and find matrix of
         % fitted values (variable Yhat)
         % Note that in this case t1 is exactly the same of t (that is
         % Yhat are evaluated in correspondence of t
        [t1,Yhat] = ode45(diffeq,t,AB0,options);
        
    end
end