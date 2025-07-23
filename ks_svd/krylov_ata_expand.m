function [V, U, alpha, beta] = krylov_ata_expand(A, V, U, c, k)

%KRYLOV_ATA_EXPAND  Expand orthonormal bases for Lanczos bidiagonalization spaces K_k(A'A, v) and K_k(AA', Av)
% function [V, U, alpha, beta] = krylov_ata_expand(A, V, U, c, k)
%
% Out: AV_k = U_k B_k and A'U_k = V_{k+1} B_{k+1}
%   B = diag(alpha) + diag(beta, 1)
%   alpha: diagonal elements of B_k
%   beta : superdiagonal elements of B_k
%
% c is for the first orthogonalization step: c = U_{m-1}'AV_m
% may be a multiple of e_m (krylov_ata) but not necessarily so (krylov_schur_svd)
%
% See also KRYLOV_ATA, KRYLOV_SCHUR_SVD
%
% Revision date: May 28, 2016
% (C) Michiel Hochstenbach 2016

if nargin < 5 || isempty(k), k = 10; end
m = size(V,2);
V = [V zeros(size(V,1), k)]; U = [U zeros(size(U,1), k)];  % Preallocation
alpha = zeros(1,k); beta  = zeros(1,k);

for j = m:(k+m-1)
  if j == m
    r = mv(A, V(:,j), 0) - U(:,1:j-1)*c;
  else
    r = mv(A, V(:,j), 0) - beta(j-m)*U(:,j-1);
  end
  r = r - U(:,1:j-1)*(U(:,1:j-1)'*r);  % Reorthogonalization
  alpha(j-m+1) = norm(r);
  U(:,j) = r / alpha(j-m+1);
  r = mv(A, U(:,j), 1) - alpha(j-m+1)*V(:,j);
  r = r - V(:,1:j)*(V(:,1:j)'*r);      % Reorthogonalization
  beta(j-m+1) = norm(r);
  V(:,j+1) = r / beta(j-m+1);
end
