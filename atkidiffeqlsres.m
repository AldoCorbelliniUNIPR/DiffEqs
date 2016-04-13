%% new version 1


%%%%%%%%%%%%%  least squares
function out=atkidiffeqlsres
% 4 thetas


% calcolare R
R= random('Normal',0,0.01,10,2);

t=(0:0.1:5)';
n=length(t);
thetaTRUE=[1 2];
v=3;
    A0=0;
etatrue=findyhat(thetaTRUE,t);
    Y=bsxfun(@plus,etatrue,0.1*randn(n,v));
theta0=[0.9 1.7];

loglik=@likfx;

[thetahat,fval,exitflag]  = fminsearch(loglik,theta0);
    
disp(thetahat)

    function dZ=likfx(theta)
    
        th2=theta(2);
        th1=theta(1);
        
        yhat=th2/(th1+th2)+(A0-th2/(th1+th2))*exp(-(th1+th2)*t);
        
        R=zeros(n,v);
        for j=1:v
            R(:,j)=Y(:,j)-yhat;
        end
        dZ=det(R'*R);
    end


  function yhat=findyhat(theta,t)
   
        th2=theta(2);
        th1=theta(1);
        
       yhat= th2/(th1+th2)+(A0-th2/(th1+th2))*exp(-(th1+th2)*t);
      
    end
out=thetahat;
end