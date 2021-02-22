clc
close all

% Reproduces Example 3.1 from Gugercin and Drmac 2015
% Generates Figure 1 from the paper

f = @(t,mu) 10*exp(-mu.*t).*(cos(4*mu.*t) + sin(4*mu.*t));
t  = linspace(1,6,10000)';
mu = linspace(0,pi,100)';

% Compute the snapshots
[T,M] = meshgrid(t,mu);
F = f(T,M);

[Uf,sf,~] = svd(F',0);

klst = [10, 20];
p = 10;
figure,
rect = [0,0, 10, 5];
set(gcf, 'Units', 'inches'); 
set(gcf, 'OuterPosition',rect);
set(gcf, 'Position', rect); 

for j = 1:length(klst)
    k = klst(j);

    % oversampling
    Uh = rangefinder(F',k,p,0);
    U = Uf(:,1:k);
    theta = subspace_angles(U,Uh);

    
    %% Compute and plot the errors
    n = 200;
    mu_eval = linspace(0,pi,n)';
    [T,M] = meshgrid(t,mu_eval);
    ft = f(T,M)';
    ftnrm = sqrt(sum(ft.^2,1));

    % DEIM Accuracy
    [pd,errd] = subsetselection(U,  'pqr');
    tp = t(pd);
    [T,M] = meshgrid(tp,mu_eval);
    fd = f(T,M)';
    f_deim = U*(U(pd,:)\fd);
    deim_err = sqrt(sum((ft-f_deim).^2,1))./ftnrm;

    % R-DEIM Accuracy
    [pr,errr] = subsetselection(Uh, 'pqr');
    tp = t(pr);
    [T,M] = meshgrid(tp,mu_eval);
    fr = f(T,M)';
    f_rdeim = Uh*(Uh(pr,:)\fr);
    rdeim_err = sqrt(sum((ft-f_rdeim).^2,1))./ftnrm;


    %% Compute and plot the bounds
    uft = U*(U'*ft); res = ft - uft;
    ftnrm  = sqrt(sum(ft.^2,1));
    Resnrm = sqrt(sum(res.^2,1));
    deim_bnd   = errd*(Resnrm)./ftnrm;
    rdeim_bnd1 = errr*(Resnrm+sin(max(theta)).*ftnrm)./ftnrm;

    I = speye(size(U,1));
    Sd = full(I(:,pd));
    Sr = full(I(:,pr));
    psi = subspace_angles(Sd,Sr);
    
    
    
    rdeim_bnd2 = (errd*(Resnrm) + errd*errr*sin(max(psi))*Resnrm ...
                + errd*errr*sin(max(theta)).*ftnrm )./ftnrm;
    
    subplot(1,2,j)
    semilogy(1:n,deim_err,'k.-',1:n,rdeim_err,'r--',1:n,deim_bnd, 'b-', 1:n, rdeim_bnd1,'r:', 1:n, rdeim_bnd2, 'm-', 'Linewidth',2)
    
    set(gca, 'FontSize',16)
    xlabel('Index','FontSize',18)
    if j == 1
        ylabel('Relative Error','FontSize',18)
        legend({'DEIM','RDEIM','DEIM - bound', 'RDEIM - bound 1', 'RDEIM - bound 2'}, 'location', 'southeast')
    end
    title(strcat('Target rank r = ',num2str(k)), 'FontSize', 20);

    sin(max(theta))
end
print -depsc figs/rdeim_error
