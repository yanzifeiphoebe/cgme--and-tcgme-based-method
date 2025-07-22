clc;      
clear all;
close all;


maxit = 2000;
%%%%%%%%%%%%%% 1-D %%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = 10000; 
% [K,b,x_true] = baart(n); 
% noise  = randn(size(b)); 
% delta=0.01;
% sigma = delta*(norm(b)/norm(noise)); 
% N = sigma*noise; 
% b = b + N;
% L = get_l(n,1);
%
% % %%%%%%%%%%%%%%%%%% 2-D  RestoreTools %%%%%%%
% N = 256;
N = 128;
K = blur(N,16,2); 
load mri 
x_true = double(D(1:N,1:N,15));
% load GaussianBlur440
% load GaussianBlur422 
% load GaussianBlur420 
% load VariantMotionBlur_small
% load satellite 
% load grain 
% load AtmosphericBlur50 
% load AtmosphericBlur30
% load AtmosphericBlur10 
% load VariantGaussianBlur1 load
% VariantGaussianBlur2 
% load VariantGaussianBlur3
% K = psfMatrix(PSF);
% x_true=im2double(f_true);
% x_true=x_true(1:N,1:N);
x_true=x_true(:);
% b=g(:);
btrue = K*x_true;
delta=0.1;
e =  randn(size(btrue(:)))/norm(randn(size(btrue(:))));
b = btrue(:) + norm(btrue(:))*delta*e;
% %  %%%%%%%%%%%%%%%% 2-D IRtools %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % [A, btrue, x_true, ProbInfo] = PRblurspeckle; [b, NoiseInfo] =
% % PRnoise(btrue, 'gauss', delta);
% %%%%%%%%%%%%%%%%%%%%%%%% L matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I = speye(N);
L1=get_l(N,1);
L = [kron(I,L1); kron(L1,I)];

%%%%%%%%%%%%%%%%%%%%%%%%% JBDQR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t1 = 0;
% tic;
% options = HyBRset('InSolv', 'none', 'x_true', x_true(:),'Iter', maxit);
% [x_out, output] = YangJLBDHyBR(K, L, b, maxit, options);
% t1=t1+toc;
% [a_jbdqr,b_jbdqr] = min(output.err);
% 
% %%% CGME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t2=0;
% tic;
% [X,V_cgme,U_c,B_c] = yangcgme(K,b,maxit,1);
% [x_cgme,err_cgme,err_freeLcg,rho_cgme,eta_cgme] = YangfunTik(K, b,x_true,X,V_cgme,U_c,B_c,L,maxit);
% t2 = t2+toc;
% [a_cgme,b_cgme] = min(err_cgme);
% 
% % %%%%%%%% TCGME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t3 = 0;
% tic
% [X_t,V_tcgme,U_t,B_t] = yangtcgme(K,b,maxit,1);
% [x_tcgme,err_tcgme,err_freeLtcg,rho_tcgme,eta_tcgme] = YangfunTik(K,b,x_true,X_t,V_tcgme,U_t,B_t,L,maxit);
% t3 = t3+toc;
% [a_tcgme,b_tcgme] = min(err_tcgme);

%%%%%%%%%%%% MTRSVD %%%%%%%%%%%%%%%%%%%%%
% generate Q and RSVD 
t4 = 0;
tic
omega = randn(N^2,maxit+6); %for 2D problems
% omega = randn(n, maxit+6); %for 1D problems
Y = K*omega;
[Q,R] = qr(Y,0);
B = (K'*Q)';
[W, S, V] = svds(B,maxit); 
U = Q*W;
[x_rm,x_s,z,f,flag,relres,iter] = YangfunMtrsvd(U,S,V,L,b,maxit);

%Compute the mtrsvd 
for i = 1:maxit
    err_rm(i) = norm(L*(x_rm(:,i)-x_true))/norm(L*x_true);
    eta_rm(i) = norm(L*x_s(:,i)-L*z(:,i));
    err_rmL(i) = norm(x_rm(:,i)-x_true)/norm(x_true);
    rho_rm(i) = norm(f(:,i)-b-K*z(:,i));
end
t4 = t4+toc;
[a_mtrsvd,b_mtrsvd] = min(err_rm);

%%%%%%%%%%%%%5 record the results %%%%%%%%%%%%%%

% time1D = {'mri' maxit delta t2 t3 t1 t4};
% xlswrite( 'C:\Users\phoebe_yang\Desktop\time1D.xlsx',time1D, 'A37:G37');
%  
% accuracy1D = {'mri' maxit delta a_jbdqr b_jbdqr a_cgme b_cgme a_tcgme b_tcgme a_mtrsvd b_mtrsvd};
% xlswrite( 'C:\Users\phoebe_yang\Desktop\accuracy1D.xlsx', accuracy1D,'A37:K37')

%%%%%%%%%%%%% only record the results of mri %%%%%%%%%%%%%%
accuracy_svds = {'mri' maxit delta a_mtrsvd b_mtrsvd};
xlswrite( 'C:\Users\phoebe_yang\Desktop\accuracy_svds.xlsx', accuracy_svds,'A14:E14')

time_svds = {'mri' maxit delta t4};
xlswrite( 'C:\Users\phoebe_yang\Desktop\time_svds.xlsx',time_svds, 'A14:D14');
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% semilogy(output.err,'.r-');hold on
% semilogy(err_cgme,'.b-');hold on
% semilogy(err_tcgme,'.g-');hold on
% semilogy(err_rm,'.m-');hold on
% legend('JBDQR','hyb-cgme','hyb-tcgme','mtrsvd')
% legend('mtrsvd')
% title('relative error'+
% axes('position',[0.4 0.2 0.5 0.4]); 
% semilogy(output.err,'.r-');hold on 
% semilogy(err_tcgme,'.g-');hold on
% % semilogy(err_rm,'.m-');hold on
% legend('jbdqr','hyb-tcgme')

% figure;
% semilogy(output.err2,'.r-');hold on
% semilogy(err_freeLcg,'.b-');hold on
% semilogy(err_freeLtcg,'.g-');hold on
% semilogy(err_rmL,'.m-');hold on
% legend('JBDQR','hyb-cgme','hyb-tcgme','mtrsvd')
% % legend('mtrsvd')
% title('relative error')

% axes('position',[0.4 0.2 0.5 0.4]); 
% semilogy(output.err2,'.r-');hold on 
% semilogy(err_freeLtcg,'.g-');hold on
% % semilogy(err_rm,'.m-');hold on
% legend('jbdqr','hyb-tcgme')

% %%%%%%%%%%%%%%%%%%%%%% L-curve of JBDQR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure;
% [reg_c,rho_c,eta_c] = l_corner(output.Rnrm,output.Xnrm);
% plot_lc(output.Rnrm,output.Xnrm,'o',2);hold on
%                
% ax = axis;hold on
% loglog([min(output.Rnrm)/100,rho_c],[eta_c,eta_c],':r',...
% [rho_c,rho_c],[min(output.Xnrm)/100,eta_c],':r')
% title(['jdbqr L-curve',' corner at ', num2str(reg_c), ' the relative error is ',num2str(output.err(reg_c))]);
% axis(ax);hold off
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% L-curve of hyb-CGME %%%%%%%%%%%%%
% figure; [reg_c,rho_c,eta_c] = l_corner(rho_cgme,eta_cgme); %
% plot_lc(rho_cgme,eta_cgme,'o',2);hold on 
% 
% ax = axis;hold on
% loglog([min(rho_cgme)/100,rho_c],[eta_c,eta_c],':r',... %
% [rho_c,rho_c],[min(eta_cgme)/100,eta_c],':r') 
% title(['hyb-cgme L-curve','corner at ', num2str(reg_c), 'the relative error is ', num2str(err_cgme(reg_c))]); axis(ax);hold off

%%%%%%%%%%%%%%%%%% L-curve of hyb-TCGME %%%%%%%%%
% figure; [reg_c,rho_c,eta_c] = l_corner(rho_tcgme,eta_tcgme);
% plot_lc(rho_tcgme,eta_tcgme,'o',2);hold on
% 
% ax = axis;hold on 
% loglog([min(rho_tcgme)/100,rho_c],[eta_c,eta_c],':r',...
% [rho_c,rho_c],[min(eta_tcgme)/100,eta_c],':r') 
% title(['hyb-tcgme L-curve',' corner at ', num2str(reg_c), 'the relative error is ', num2str(err_tcgme(reg_c))]); axis(ax);hold off
% %
% %%%%%%%%%%%%%%%% L-curve of MTRSVD %%%%%%%%%%
% figure; [reg_c,rho_c,eta_c] = l_corner(rho_rm,eta_rm); 
% plot_lc(rho_rm,eta_rm,'o',2);hold on 
% ax = axis;hold on 
% loglog([min(rho_rm)/100,rho_c],[eta_c,eta_c],':r',... 
% [rho_c,rho_c],[min(eta_rm)/100,eta_c],':r') 
% title(['mtrsvd L-curve','corner at ', num2str(reg_c), ' the relative error is ', num2str(err_rm(reg_c))]); axis(ax);hold off
% 
% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %for 1D problems
% figure; 
% plot(x_true,'.r-'); hold on
% plot(x_cgme(:,b_cgme),'.m-');hold on
% plot(x_tcgme(:,b_tcgme),'.c-');hold on
% plot(x_out(:,b_jbdqr),'.g-');hold on
% plot(x_rm(:,b_mtrsvd),'.b-');hold on
% legend('original image','hyb-cgme','hyb-tcgme','JBDQR','MTRSVD') 

% % for 2D problems
% % figure; 
% subplot(2,3,1) 
% imagesc(reshape(x_true,N,N)) 
% colormap gray, axis image off 
% title('original image') 
% subplot(2,3,2)
% imagesc(reshape(x_cgme(:,b_cgme),N,N)) 
% colormap gray, axis image off
% title('hyb-cgme') 
% subplot(2,3,3) 
% imagesc(reshape(x_tcgme(:,b_tcgme),N,N))
% colormap gray, axis image off 
% title('hyb-tcgme') 
% subplot(2,3,4)
% imagesc(reshape(x_out(:,b_jbdqr),N,N)) 
% colormap gray, axis image off
% title('JBDQR')
% subplot(2,3,5)
% imagesc(reshape(x_rm(:,b_mtrsvd),N,N)) 
% colormap gray, axis image off
% title('MTRSVD')
% subplot(2,3,6) 
% imagesc(reshape(b,N,N)) 
% colormap gray, axis image off 
% title('blur image') 
% 
% 
%for only mtrsvd algorithm
figure;subplot(1,3,1)
imagesc(reshape(x_rm(:,b_mtrsvd),N,N)) 
colormap gray, axis image off
title('MTRSVD')
subplot(1,3,2)
imagesc(reshape(b,N,N)) 
colormap gray, axis image off
title('blurred image')
subplot(1,3,3) 
imagesc(reshape(x_true,N,N)) 
colormap gray, axis image off 
title('original image') 
