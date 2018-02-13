function [x,y,w,z,u,it,error_abs,error_rel] = ...
    multi_dimensional_HITS(A,x0,w0,u0,alpha,b,tol)
% Implementation of the iterative method for the computation of the 
% multi-dimensional HITS centrality vectors introduced in the paper
% =========================================================================
% Multi-dimensional HITS:                                                 
% An always computable ranking for temporal multi-layer directed networks
% by Francesca Arrigo and Francesco Tudisco
% =========================================================================
% PLEASE CITE OUR PAPER IF YOU USE THIS CODE
%
%
%
% This implementation of the method requires the Matlab Tensor Toolbox:
% @misc{TTB_Software,
%  author = {Brett W. Bader and Tamara G. Kolda and others},
%  title = {MATLAB Tensor Toolbox Version 2.6},
%  howpublished = {Available online},
%  month = {February},
%  year = {2015},
%  url = {http://www.sandia.gov/~tgkolda/TensorToolbox/}
% }
%
%
%
% INPUT
% =========================================================================
% A     :   Adjacency tensor (sparse) in the sptensor format of dimension
%           nV x nV x nL x nL x NT
% x0    :   Starting vector of size nV
% w0    :   Starting vector of size nL
% u0    :   Starting vector of size nT
% alpha :   Vector of exponents
% b     :   Eigenvector of the martix M defined in terms of alpha as in the
%           paper. Recall that b is required to be such that sum(b) = 1
% tol   :   Tolerance for the stopping criterion
%
% OUTPUT
% =========================================================================
% x     :   Hub centrality of nodes
% y     :   Authority centrality of nodes
% w     :   Broadcasting centrality of layers
% z     :   Receiving centrality of layers
% u     :   Time centrality
% error_abs   :  Absolute error between each two consecutive iterations with
%                respect to the norm ||x_k - x_{k+1}||_{b, infty} defined
%                in the paper
% error_rel   :  Relative error between each two consecutive iterations  
%                i.e. ||x_k - x_{k+1}||_{b,infty}/||x_{k+1}||_{b,infty}
%
%
%
% This version is dated: 10 Feb. 2018  [FA]


%%% Check that rho(M) < 1
M = ones(5,1)*alpha - diag(alpha);
[~,D] = eig(M); D = diag(D); 
[rho,~] = max(D);
if rho >=1 
    warning('rho >=1, the method might not converge!');
end


x = x0; %./norm(x0,'inf');
y = x0;
w = w0; %/norm(w0,1);
z = w0;
u = u0;


error_rel = 1; error_abs = 1;
it = 0;


while error_rel(end) > tol
    xold = x;
    yold = y;
    wold = w;
    zold = z;
    uold = u;
    

    xx = ttv(A,{y,w,z,u},[2 3 4 5]); 
    xx = double(tenmat(xx,1)).^alpha(1);
    %%% The following might be faster:
    % xx = ttv(A,{y,w,z},[2 3 4]);
    % xx = double(sptenmat(xx,1))*u;
    % xx = xx.^alpha(1);
    
    
    
    yy = ttv(A,{x,w,z,u},[1 3 4 5]);
    yy = double(tenmat(yy,1)).^alpha(2);
    %%% The following might be faster:
    % yy = ttv(A,{x,w,z},[1 3 4]);
    % yy = double(sptenmat(yy,1))*u;
    % yy = yy.^alpha(2);

    
    ww = ttv(A,{x,y,z,u},[1 2 4 5]);
    ww = double(tenmat(ww,1)).^alpha(3);
    %%% The following might be faster:
    % ww = ttv(A,{x,y,z},[1 2 4]);
    % ww = double(sptenmat(ww,1))*u;
    % ww = ww.^alpha(3);
    
    zz = ttv(A,{x,y,w,u},[1 2 3 5]);
    zz = double(tenmat(zz,1)).^alpha(4);
    %%% The following might be faster:
    % zz = ttv(A,{x,y,w},[1 2 3]);
    % zz = double(sptenmat(zz,1))*u;
    % zz = zz.^alpha(4);
    
    uu = ttv(A,{x,y,w,z},[1 2 3 4]);
    uu = double(tenmat(uu,1)).^alpha(5);
    %%% The following might be faster: 
    %uu = ttv(A,{x,y,w},[1 2 3]);
    % uu = z'*double(sptenmat(uu,1));
    % uu = uu(:);
    % uu = uu.^alpha(5);
    
    x = abs(xx)/norm(abs(xx),'inf'); 
    y = abs(yy)/norm(abs(yy),'inf');
    w = abs(ww)/norm(abs(ww),'inf');
    z = abs(zz)/norm(abs(zz),'inf');
    u = abs(uu)/norm(abs(uu),'inf');
 
    
    new_err = b(1)*norm(xold-x,'inf') + b(2)*norm(yold-y,'inf') + ...
        b(3)*norm(wold-w,'inf') + b(4)*norm(zold-z,'inf') + ...
        b(5)*norm(uold-u,'inf');
    error_abs = [error_abs;new_err];
    new_err = new_err./ (b(1)*norm(x,'inf') + b(2)*norm(y,'inf') + ...
        b(3)*norm(w,'inf') + b(4)*norm(z,'inf') + ...
        b(5)*norm(u,'inf'));
    error_rel = [error_rel; new_err];
    

    
    it = it + 1;
    
end     
