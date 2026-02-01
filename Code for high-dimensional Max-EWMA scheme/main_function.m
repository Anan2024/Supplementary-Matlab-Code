% Algorithm 1: Searching for the UCL value to achieve target MRL
% Input parameters
MRL_0 = 250;       % Desired in-control median run length
r = 50000;          % Number of simulation replications
T_0 = 800;        % Maximum run length to simulate (T_0 > MRL_0)
m = 100;            % Size of reference sample
n = 5;             % Size of test sample
lambda = 0.05;      % Smoothing parameter for EWMA control chart
d = 2;             % Distance parameter for control limits
p=3;               % dimensionality


tic;
UCL = bisection_method_to_refine_UCL(MRL_0, r, T_0, m, n, lambda, d,p);
elapsed_time = toc;



    
