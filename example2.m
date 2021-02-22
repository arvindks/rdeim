clc
close all

% This code generates figure 5 from the paper


%% 2D function
g1 = @(x,y,mu) 1./sqrt( ((1-x)-(0.99*mu(1)-1)).^2 + ...
                    ((1-y)-(0.99*mu(2)-1)).^2 + 0.1.^2);
                
s = @(x,y,mu) g1(x,y,mu) + g1(1-x,1-y,1-mu) + ...
                g1(1-x,y,[1-mu(1),mu(2)]) + g1(x,1-y,[mu(1),1-mu(2)]);

%Domain 
nx = 25;    ny = 25;
x = linspace(0.0,1.0,nx);   y = linspace(0.0,1.0,ny);
[X,Y] = meshgrid(x,y);
pts = [X(:),Y(:)];close

%Parameter
mu1 = linspace(0,1,25);
mu2 = linspace(0,1,25);
[Mu1,Mu2] = meshgrid(mu1,mu2);
mu = [Mu1(:),Mu2(:)];

% Collect snapshots
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


%% POD-DEIM approximation for each dimension size
ns = 30;    %Number of snapshots

Errnrm1 = zeros(ns,1);  Errnrm2 = zeros(ns,1);  Errnrm3 = zeros(ns,1);  Errnrm4 = zeros(ns,1);
Errcst1  = zeros(ns,1); Errcst2  = zeros(ns,1); Errcst3  = zeros(ns,1); Errcst4  = zeros(ns,1);


pqrtime = zeros(ns,1);
hybridtime = zeros(ns,1);
lstime = zeros(ns,1);

Fnrm = max(sqrt(sum(F.^2,1)));
[u,~,~] = svd(S,0);


% Number of repetitions for timing
nreps = 10;

for i = 2:ns
    ur = u(:,1:i);
    
    % PQR
    tic
    for j = 1:nreps
        [p,err] = subsetselection(ur, 'pqr');
    end
    pqrtime(i) = toc/nreps; 
    Errcst1(i) = err;
    Fh = ur*(ur(p,:)\F(p,:));
    
    Err = F-Fh;
    Errnrm1(i) = max(sqrt(sum(Err.^2,1)))./Fnrm;
    
    % RRQR
    tic
    for j = 1:nreps
        [pr,~,err] = hybrid(ur,'det','pqr','method','const');
    end
    hybridtime(i) = toc/nreps;  
    Errcst2(i) = err;
    
    Fh = ur*(ur(pr,:)\F(pr,:));
    
    
    % LS sampling
    
    Err = F-Fh;
    Errnrm2(i) = max(sqrt(sum(Err.^2,1)))./Fnrm; 
 
    tic
    for j = 1:nreps
        [pr,d,err] = randls(ur,'method','const');
    end
    lstime(i) = toc/nreps;
    
    Errcst4(i) = err;   s = length(d);  D = spdiags(d,1:1,s,s);
    Fh = ur*(D*ur(pr,:)\(D*F(pr,:)));
    
    Err = F-Fh;
    Errnrm4(i) = max(sqrt(sum(Err.^2,1)))./Fnrm; 
    
end

%% 
figure, 
semilogy(2:ns,Errnrm1(2:ns),2:ns,Errnrm2(2:ns),'--',2:ns,Errnrm4(2:ns),':','LineWidth',4.0)
set(gca, 'FontSize',18)
legend('PQR','Hybrid',  'LS') % 'LS \epsilon = 0.99, \delta = 0.1',
xlabel('Target rank','FontSize',18)
title('Relative Error','FontSize',20)
print 'figs/example2_err' -depsc


figure, semilogy(2:ns, pqrtime(2:ns), ...
    2:ns, hybridtime(2:ns), '--', 2:ns, lstime(2:ns), ':','LineWidth',4.0)
set(gca, 'FontSize',18)
legend({'PQR','Hybrid',  'LS'}, 'location', 'southeast') % 'LS \epsilon = 0.99, \delta = 0.1',
xlabel('Target rank','FontSize',18)
title('Computational Time','FontSize',20)
print 'figs/example2_time' -depsc

figure, 
semilogy(2:ns,Errcst1(2:ns),2:ns,Errcst2(2:ns),'--',2:ns,Errcst4(2:ns),':','LineWidth',4.0)
set(gca, 'FontSize',18)
legend({'PQR','Hybrid', 'LS'}, 'location','best'); %'LS \epsilon = 0.99, \delta = 0.1',
xlabel('Target rank','FontSize',18)
title('Error constant','FontSize',20)
print 'figs/example2_errc' -depsc

figure,
[pr,pls,err] = hybrid(ur,'det','rrqr','method','const');
plot(mu(pls,1),mu(pls,2),'rx','MarkerSize',7.0), hold on
plot(mu(pr,1),mu(pr,2),'ko','MarkerSize',8.0)
xlabel('x_1'); ylabel('x_2')
set(gca, 'FontSize',18)
legend('LS','Hybrid')
title('Hybrid vs LS','FontSize',20)
print 'figs/example2_hybridls' -depsc

figure,
plot(mu(p,1),mu(p,2),'bs','MarkerSize',7.0,'MarkerFaceColor','b'), hold on
plot(mu(pr,1),mu(pr,2),'ro','MarkerSize',8.0)
xlabel('x_1'); ylabel('x_2')
set(gca, 'FontSize',18)
legend('sRRQR','Hybrid')
title('Hybrid vs sRRQR', 'FontSize',18)
print 'figs/example2_hybridqr' -depsc


