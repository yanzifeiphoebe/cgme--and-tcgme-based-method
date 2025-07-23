function r = rand1(n, p)

%RAND1  Creates a random vector with norm = 1
% function r = rand1(n, p)
%
% See also RANDN, RANDORTH
%
% Revision date: April 2, 2020
% (C) Michiel Hochstenbach 2022

if nargin < 2 || isempty(p), p = 2; end

if p == 2
  r = randn(n,1);  r = r / norm(r);  % This is randomly distributed on S^(n-1)
elseif p == 1
  r = randn(n,1);  r = r / norm(r, 1);
else % p = Inf
  r = randn(n,1);  r = r / norm(r, Inf);
end
