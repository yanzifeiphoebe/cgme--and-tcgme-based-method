function [X,V,U,B] = yangtcgme(A,b,k,reorth)
%MCGME_B Solution of least squares problems by Lanczos bidiagonalization.
%
% [X,rho,eta] = yangmcgme_b(A,b,k,reorth,s)
%
% Performs k steps of the MCGME Lanczos bidiagonalization algorithm
% applied to the system
%    min || A x - b || .
% The routine returns all k solutions, stored as columns of
% the matrix X.
%
% The solution norm are returned in eta
% The residual norm are returned in rho
%
%
%
% Reorthogonalization is controlled by means of reorth:
%    reorth = 0 : no reorthogonalization (default),
%    reorth = 1 : reorthogonalization by means of MGS.

% Reference: C. C. Paige & M. A. Saunders, "LSQR: an algorithm for
% sparse linear equations and sparse least squares", ACM Trans.
% Math. Software 8 (1982), 43-71.

%

% The fudge threshold is used to prevent filter factors from exploding.


% Initialization.
if (k < 1), error('Number of steps k must be positive'), end
if (nargin==3), reorth = 0; end
%if (nargout==4 && nargin<5), error('Too few input arguments'), end
[m,n] = size(A); X = zeros(n,k);
if (reorth==0)
    UV = 0;
elseif (reorth==1)
    U = zeros(m,k); V = zeros(n,k); UV = 1;
    if (k>=n), error('No. of iterations must satisfy k < n'), end
else
    error('Illegal reorth')
end
if (nargout > 1)
    eta = zeros(k,1); rho = eta;
    
end


% Prepare for LSQR iteration.
v = zeros(n,1); beta = norm(b);
if (beta==0), error('Right-hand side must be nonzero'), end
u = b/beta; if (UV), U(:,1) = u; end
r = A'*u; alpha = norm(r);   % A'*u;
v = r/alpha; if (UV), V(:,1) = v; end
B(1,1) = alpha; beta_old = beta;

% zeta = -1;

% Perform Lanczos bidiagonalization with/without reorthogonalization.
for i=1:k+1
    
    % Update the solution.
%     zeta = (-beta/alpha)*zeta;
%     x = x + zeta*v; X(:, i) = x;
%     eta(i) = norm(x);
    
    %%%%%%%%%%%%%%  Compute A*v - alpha*u.  %%%%%%%%%%%%%%%%%%
    p = A*v - alpha*u;
    if (reorth==0)
        beta = norm(p); u = p/beta;
    else
        for j=1:i, p = p - (U(:,j)'*p)*U(:,j); end
        beta = norm(p); u = p/beta;
    end
    B(i+1,i) = beta;
    
    %%%%%%%%%%%%%%  Compute A'*u - beta*v.  %%%%%%%%%%%%%%%%
    r = A'*u - beta*v;
    if (reorth==0)
        alpha = norm(r); v = r/alpha;
    else
        for j=1:i, r = r - (V(:,j)'*r)*V(:,j); end
        alpha = norm(r); v = r/alpha;
    end
    B(i+1,i+1) = alpha;
    
    %%%%%%%%% Store U and V if necessary.%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (UV), U(:,i+1) = u; V(:,i+1) = v; end
    
    %%%%%%%% Computing the SVD and the best rank-k approximation of B. %%%%
    [P,S,Q] = svd(B); 
    y = beta_old*Q(:,1:i)*(diag(diag(S(1:i,1:i)))\P(1,1:i)');
    x = V(:,1:i+1)*y; X(:,i) = x; eta(i) = norm(x);
    r = b - A*x; rho(i) = norm(r);
end