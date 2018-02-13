%%%%%%%%%%%%%%%%
%%% This file is a guideline example to the use of multi-dimensional HITS
%%%%%%%%%%%%%%%%

%%% The current implementation of Algorithm 1 for the computation of  
%%% the multi-dimensional HITS centrality requires the Tensor Toolbox:
%%% @misc{TTB_Software,
%%%  author = {Brett W. Bader and Tamara G. Kolda and others},
%%%  title = {MATLAB Tensor Toolbox Version 2.6},
%%%  howpublished = {Available online},
%%%  month = {February},
%%%  year = {2015},
%%%  url = {http://www.sandia.gov/~tgkolda/TensorToolbox/}
%%% }

%%% Adjacency tensor of the network. 
%%% The tensor used in the paper (temporal citation network) can be parsed
%%% from the txt files from the folder citation_network_data

%%% This guide uses a sparse binary random tensor
nV = 30; % number of nodes
nL = 15;  % number of layers
nT = 10;  % number of time stamps

T = sptenrand([nV nV nL nL nT], .05); T = sptensor(T > 0); % binary random tensor

%%% Model parameters alpha, matrix M and associated Perron eigenvector b
%%% Here we use alpha_1 = ... = alpha_5 = 1/5
%%% Different choices of alpha require to check the condition rho < 1
%%% A random choice of alpha and b can be obtained with the function
%%% [alpha,b] = rand_alpha(1);
alpha = ones(5,1)./5; 
M = ones(5,1)*alpha - diag(alpha);
[V,D] = eig(M); D = diag(D); 
[rho,ind] = max(D); % check rho < 1 if other alpha are used
b = abs(V(:,ind));
b = b./norm(b,1);

tol = 1e-06; % tolerance

%%% Starting vectors. Here we use the all-ones vectors
x0 = ones(nV,1); 
w0 = ones(nL,1);
u0 = ones(nT,1);

[h,a,b,r,t,...
 num_iterations,...
 error,...
 relative_error] = multi_dimensional_HITS(T,x0,w0,u0,alpha,b,tol);

                                        
                                        
