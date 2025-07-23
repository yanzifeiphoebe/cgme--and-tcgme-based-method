function e = element(A, i, j, n)

%ELEMENT  Element A(i,j) of matrix or function
% function e = element(A, i, j, n)
%
% Practical to avoid incorrect syntax such as (A*B)(2,1), or when A is a function
% n is only needed when A is a function A: C^n -> C^m
%
% Revision date: June 10, 2020
% (C) Michiel Hochstenbach 2020

if nargin < 2, i = []; end
if nargin < 3, j = []; end

if isnumeric(A)        % Matrix
  if min(size(A)) > 1
    if isempty(i), e = A(:,j);
    elseif isempty(j), e = A(i,:); else e = A(i,j); end
  else, e = A(i); end  % Vector
else                   % Function
  if isempty(j), error('j has to be nonempty when A is a function'); end
  e = mv(A,unv(j,n));
  if ~isempty(i), e = e(i); end
end
