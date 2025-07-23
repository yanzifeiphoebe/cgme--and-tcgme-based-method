function [V, U, alpha, beta] = krylov_ata(A, v1, k, full, reortho)

%KRYLOV_ATA  Construct orthonormal bases for Lanczos bidiagonalization spaces K_k(A'A,v) and K_k(AA',Av)
% function [V, U, alpha, beta] = krylov_ata(A, v1, k, full, reortho)
%   
% Out: AV_k = U_k B_k and A'U_k = V_{k+1} B_{k,k+1}'
% alpha    diagonal elements of bidiagonal matrix B_k
% beta     superdiagonal elements of B_k
%
% k        dimension of Krylov spaces
% full     0: basis for V_k,     length(beta) = k
%          1: basis for V_{k+1}, length(beta) = k+1
% reortho  0: no reorthogonalization (not recommended)
%          1: reorthogonalization of V vectors
%          2: reorthogonalization of U vectors and V vectors (default)
%
% Reference (for instance) Golub & Van Loan, 4th ed., 2013
%
% See also KRYLOV_ATA_EXPAND, KRYLOV_SYM, KRYLOV_SCHUR_SVD
%
% Revision date: May 13, 2022
% (C) Michiel Hochstenbach 2022

if nargin < 2 || isempty(v1)
  if ~isnumeric(A)
    error('Need a starting vector if A is a function handle')
  end
  v1 = rand1(size(A,2));
end
if nargin < 3 || isempty(k),       k = 10; end
if nargin < 4 || isempty(full),    full = 1; end
if nargin < 5 || isempty(reortho), reortho = 2; end

alpha = zeros(1,k); beta  = zeros(1,k-1); % Preallocation
if reortho
  V = zeros(length(v1), k); V(:,1) = v1 / norm(v1);
else
  v = v1 / norm(v1);
end

for j = 1:k
  if reortho
    r = mv(A, V(:,j), 0);
    if j == 1 && reortho == 2, U = zeros(length(r), k); end % Preallocation
  else
    r = mv(A, v, 0);
  end
  if j > 1
    if reortho == 2
      r = r - beta(j-1)*U(:,j-1);      % Reorthogonalization of U vectors
      r = r - U(:,1:j-1)*(U(:,1:j-1)'*r);
    else
      r = r - beta(j-1)*u;
    end
  end
  alpha(j) = norm(r);
  if reortho == 2
    U(:,j) = r / alpha(j); r = mv(A, U(:,j), 1);
  else
    u = r / alpha(j); r = mv(A, u, 1);
  end
  if reortho  
    r = r - alpha(j)*V(:,j);           % Reorthogonalization of V vectors
    r = r - V(:,1:j)*(V(:,1:j)'*r);
  else
    r = r - alpha(j)*v;
  end
  if j < k || full
    beta(j) = norm(r);
    if reortho
      V(:,j+1) = r / beta(j);
    else
      v = r / beta(j);
    end
  end
end

if ~reortho,    V = v; end
if reortho < 2, U = u; end
