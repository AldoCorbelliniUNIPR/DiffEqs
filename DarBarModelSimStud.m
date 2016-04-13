%% Simulation study which enables to analyze the stability of parameter estimates
lambda=[1 1];

% N = number of trials
N=8;
 N=200;


% v=size of the data (v corresponds to m in Pascal's message)
% In this simulation study v=3 (there are three responses)
v=3;

% sigma= sigma of normal random errors to add to fitted values
% Length of sigma is equal to v because the responses can be affected by
% errors of different magnitudes
sigma=[1 1 1];
% sigma=sigma*1e-6;

% tspan = time range
tspan = [0 50];

% Define initial conditions
AB0=[1 0 0];

% theta= true parameter vector
theta=[0.7 0.2];

% theta0 = starting values for the parameters
theta0=[0.7 0.22];

% Define the tolerance
%AbsTol=[1e-10 1e-10 1e-10];
AbsTol=[1e-7 1e-7 1e-7];
options = odeset('RelTol',1e-7,'AbsTol',AbsTol);

nsimul=100;
Par=zeros(nsimul,2);


for jj=1:nsimul
    % d = vector with 8 numbers which are 4 or 5 or 6 with the same probability
    d=ceil(rand(N,1)*3+3);
    
    Tall=NaN(max(d),N);
    for j=1:N
        Tall(1:d(j),j)=tspan(2)*rand(d(j),1);
    end
    
    % Obtain simulated data
    % Each column refers to a particular simulated trial
    
    % Ysim is a 3D array (the third dimension is associated with v = number of
    % responses)
    Ysim=zeros(max(d),N,v);
    
    
    % etaA, etaB and etaC will be the expected response function in
    % correspondence of times Tall. To the expected responses we add a normal
    % random noise whose variance is sigma^2
    j=1;
    etaA=AB0(1)*exp(-theta(1)*Tall);
    Ysim(:,:,j)=etaA+sigma(j)*randn(max(d),N);
    
    
    j=2;
    etaB=AB0(1)*theta(1)*(exp(-theta(2)*Tall)-exp(-theta(1)*Tall))/(theta(1)-theta(2));
    Ysim(:,:,j)=etaB+sigma(j)*randn(max(d),N);
    
    j=3;
    etaC=AB0(1)-etaA-etaB;
    Ysim(:,:,j)=etaC+sigma(j)*randn(max(d),N);
    
    % SelTrial = vectors made up of integers which enables us to select the
    % trials we want analyze (if SelTrial=1:N data from all trials are taken)
    SelTrial=1:N;
    
    Y=reshape(Ysim(:,SelTrial,:),max(d)*N,3);
    t=Tall(:);
    boo=~isnan(t);
    
    % Select observations corresponding to the trials
    Y=Y(boo,:);
    t=t(boo);
    % Nobs = number of observations
    Nobs=length(t);
    
    % Note that the elements of t must be ordered
    [~,sortindexes]=sort(t);
    tsor=t(sortindexes);
    Ysor=Y(sortindexes,:);
    
    
    
    out=NLSminODE(Ysor,tsor,theta0,AB0,lambda,options);
    Par(jj,:)=sin(out);
end

%% Plottting part

subplot(1,2,1)
boxplot(Par(:,1))
title('Distribution of $\hat \theta_1$','Interpreter','Latex')

subplot(1,2,2)
boxplot(Par(:,2))
title('Distribution of $\hat \theta_2$','Interpreter','Latex')


