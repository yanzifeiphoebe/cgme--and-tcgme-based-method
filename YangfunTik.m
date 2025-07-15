
function [x_k,err_k,err_freeL,rho,eta] = YangfunTik(A,b,xtrue,X,V,U,B,L,k)

 xnorm = norm(L*xtrue);
 N=size(A,2);
% [x_k,err_k,res_k,sol_k,flag,relres,iter] = YangfunTik(A,b,xtrue,X,V,L,k)
%%%%%% Compute the lsqrTik and the relative errors %%%%%%%%%%%%%%%%%%
for i=1:k
    Q = V(:,1:i);
    d=L*X(:,i);
    z = lsqr(@(z,tflag)afun(z,L,Q,tflag),d,1e-6,N);
    x_k(:,i) = X(:,i) - z;
    err_k(i) = norm(L*(x_k(:,i)-xtrue))/xnorm;
    err_freeL(i) = norm((x_k(:,i)-xtrue))/norm(xtrue);
    rho(i,1) = norm(A*x_k(:,i)-b);
%     U_hat = U(:,1:i)'*x_k(:,i);
%     rho(i,1) = norm(B(1:i+1,1:i)*U_hat-norm(b)*eye(i+1,1));
    eta(i,1) = norm(L*x_k(:,i));
    
    
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