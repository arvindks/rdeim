clc
close all


% This code generates figures 2 and 3 from the paper 

%% 2D function
g1 = @(x,y,mu) 1./sqrt( ((1-x)-(0.99*mu(1)-1)).^2 + ...
                    ((1-y)-(0.99*mu(2)-1)).^2 + 0.1.^2);
                
s = @(x,y,mu) g1(x,y,mu) + g1(1-x,1-y,1-mu) + ...
                g1(1-x,y,[1-mu(1),mu(2)]) + g1(x,1-y,[mu(1),1-mu(2)]);

%Domain 
nx = 100;    ny = 100;
x = linspace(0.0,1.0,nx);   y = linspace(0.0,1.0,ny);
[X,Y] = meshgrid(x,y);
pts = [X(:),Y(:)];close

%Parameter
mu1 = linspace(0,1,25);
mu2 = linspace(0,1,25);
[Mu1,Mu2] = meshgrid(mu1,mu2);
mu = [Mu1(:),Mu2(:)];

%Collect snapshots
S = zeros(nx*ny,size(mu,1));
for i = 1:size(mu,1)
   sv = s(X,Y,mu(i,:));
   S(:,i) =  sv(:);
end


%% Create subdomain to test the error 
mux = linspace(0,1,11);  muy = linspace(0,1,11);
[XX,YY] = meshgrid(mux,muy);
Mutest = [XX(:), YY(:)];
F = zeros( size(pts,1), size(Mutest,1));
for i = 1:size(Mutest,1)    
    sv = s(X,Y,Mutest(i,:));
    F(:,i) = sv(:);
end



%% Adaptive range - frobenius norm 
[u,s,~] = svd(S,0);
s = diag(s);
tol = 10.^(0:-1:-6);
basisdim1 = 0*tol;
basisdim2 = 0*tol;
for i = 1:length(tol)
        [Q,~,~] = randQB_EI_auto(S, tol(i), 10, 1);
        basisdim1(i) = size(Q,2);
        [Q,errest,iter] = adaptrange(S, 40, 10, tol(i),'fro');
        basisdim2(i) = size(Q,2);
end

tolrank = 0*tol;
for i = 1:length(tol)
    tolrank(i) = find(tol(i) >= sqrt(1-cumsum(s.^2)/sum(s.^2)),1,'first');
end


figure, semilogx(tol, tolrank, tol, basisdim1, 'LineWidth', 4.0)
l = legend('r_\epsilon','Adaptive'); %, 'Adaptive - One Pass');
set(l, 'FontSize',16)
xlabel('Relative error \epsilon','FontSize',18)
ylabel('Dimension', 'FontSize',18)
title('Adaptive range finder')
set(gca, 'FontSize',18)
axis tight
print 'figs/example1_adapt_fro' -depsc


%% Testing accuracy of subspace iterations
c = 0;
errsub = zeros(length(20:5:35),5);
for k = 20:5:35
    c = c+1;
    Omega = randn(size(S,2),k+10);
    for q = 0:4
       [Q] = subspaceiteration(S,Omega,q);
       [uh,~,~] = randsvd(S,Q,k);
       theta = subspace(u(:,1:k),uh);
       errsub(c,q+1) = sin(max(theta));
    end
    
end

%
figure, semilogy(0:4,errsub','LineWidth',4.0)
legend('r=20','r=25','r=30','r=35')
xlabel('Iteration number')
ylabel('$\sin\theta_{max}$','Interpreter','LaTeX')
title('Effect of subspace iteration with p = 10')
set(gca,'FontSize',18)
print 'figs/acc_p_10' -depsc


% p = 20
c = 0;
errsub = zeros(length(20:5:35),5);
for k = 20:5:35
    c = c+1;
    Omega = randn(size(S,2),k+20);
    for q = 0:4
       [Q] = subspaceiteration(S,Omega,q);
       [uh,~,~] = randsvd(S,Q,k);
       theta = subspace(u(:,1:k),uh);
       errsub(c,q+1) = sin(max(theta));
    end
    
end

%
figure, semilogy(0:4,errsub','LineWidth',4.0)
legend('r=20','r=25','r=30','r=35')
xlabel('Iteration number')
title('Effect of subspace iteration with p = 20')
ylabel('$\sin\theta_{max}$','Interpreter','LaTeX')
set(gca,'FontSize',18)
print 'figs/acc_p_20' -depsc



%% Testing influence of oversampling

c = 0;
errsub = zeros(length(20:5:35),5);
Omega = randn(size(S,2),35+25);
for k = 20:5:35
    c = c+1;
    d = 0;
    for p = 5:5:25
       d = d+1;
       
       [Q] = subspaceiteration(S,Omega(:,1:k+p),0);
       [uh,~,~] = randsvd(S,Q,k);
       theta = subspace(u(:,1:k),uh);
       errsub(c,d) = sin(max(theta));
    end
    
end

%
figure, semilogy(5:5:25,errsub','LineWidth',4.0)
legend('r=20','r=25','r=30','r=35')
xlabel('Oversampling parameter')
ylabel('$\sin\theta_{max}$','Interpreter','LaTeX')
title('Effect of oversampling parameter with q = 0')
set(gca,'FontSize',18)
print 'figs/over_q_0' -depsc

c = 0;
Omega = randn(size(S,2),35+25);
errsub = zeros(length(20:5:35),5);
for k = 20:5:35
    c = c+1;
    d = 0;
    for p = 5:5:25
       d = d+1;
      
       [Q] = subspaceiteration(S,Omega(:,1:k+p),2);
       [uh,~,~] = randsvd(S,Q,k);
       theta = subspace(u(:,1:k),uh);
       errsub(c,d) = sin(max(theta));
    end
    
end

%
figure, semilogy(5:5:25,errsub','LineWidth',4.0)
legend('r=20','r=25','r=30','r=35')
xlabel('Oversampling parameter')
ylabel('$\sin\theta_{\max}$','Interpreter','LaTeX')
title('Effect of oversampling parameter with q = 2')
set(gca,'FontSize',18)
print 'figs/over_q_2' -depsc




