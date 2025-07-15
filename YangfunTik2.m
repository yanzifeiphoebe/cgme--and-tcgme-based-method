
function [err_k1,err_k2,err_k3] = YangfunTik2(A,b,xtrue,X,V,L,k)

 xnorm = norm(L*xtrue);
 N=size(A,2);
% [x_k,err_k,res_k,sol_k,flag,relres,iter] = YangfunTik(A,b,xtrue,X,V,L,k)
%%%%%% Compute the lsqrTik and the relative errors %%%%%%%%%%%%%%%%%%
for i=1:k
    Q = V(:,1:i);
    d=L*X(:,i);
    z1 = lsqr(@(z,tflag)afun(z,L,Q,tflag),d,1e-6,N);
    z2 = lsqr(@(z,tflag)afun(z,L,Q,tflag),d,1e-5,N);
    z3 = lsqr(@(z,tflag)afun(z,L,Q,tflag),d,1e-4,N);
    x_k1(:,i) = X(:,i) - z1;
    x_k2(:,i) = X(:,i) - z2;
    x_k3(:,i) = X(:,i) - z3;
    err_k1(i) = norm(L*(x_k1(:,i)-xtrue))/xnorm;
    err_k2(i) = norm(L*(x_k2(:,i)-xtrue))/xnorm;
    err_k3(i) = norm(L*(x_k3(:,i)-xtrue))/xnorm;
%     rho(i,1) = norm(A*x_k(:,i)-b);
%     eta(i,1) = norm(L*x_k(:,i));
    
    
end
end
function y =afun(z,A,B,transp_flag)
if strcmp(transp_flag,'transp')   % y = (A(I_n-BB^T))' * z;
    s = A'*z;
    t=B'*s;
    t=B*t;
    y = s-t;
elseif strcmp(transp_flag,'notransp') % y = (A(I_n-BB^T)) * z;
    s = A*z;
    t=B'*z;
    t=B*t;
    t = A*t;
    y = s-t;
    
end
end