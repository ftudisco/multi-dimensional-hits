function [alpha,b] = rand_alpha(k)
% Creates k 5-tuple of values such that \rho(alpha) < 1, where \rho(\alpha)
% is the leading eigenvalue of the matrix
%               M = ones(5,1)*alpha(:)^t - diag(alpha)
% The function also returns the normalized eigenvectors (in norm 1)
% associated to \rho(alpha).
% The outputs are both of size 5 x k. 

c = 0;
alpha = zeros(5,k);
b = alpha;

while c<k
    a = rand(1,5); 
    M = ones(5,1)*a - diag(a);
    [V,D] = eig(M); D = diag(D); 
    [rho,ind] = max(D);
    if rho < 1
        c = c+1;
        alpha(:,c) = a(:);
        b(:,c) = abs(V(:,ind));
        b(:,c) = b(:,c)/norm(b(:,c),1);
    end
end

