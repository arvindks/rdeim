clc
close all

% This code generates Figure 4 from the paper


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


%% POD-DEIM approximation for each dimension size
ns = 30;    %Number of snapshots

Errnrm = zeros(ns,1);
Errnrmr = zeros(ns,1);

Fnrm = max(sqrt(sum(F.^2,1)));
ratio_intersect = zeros(ns,1);

nreps = 100;

tic
for i = 1:nreps
    [u,s,~] = svd(S,0);
end
svdtime = (toc/nreps)*ones(ns,1); 



rtime = zeros(ns,1);
Dnrm = zeros(ns,1);
Dhnrm = zeros(ns,1);
for i = 1:ns
    
    
    svdtime_ = 0;
    rtime_ = 0;
    count_intersect = 0;
    for j = 1:nreps
        
        % DEIM Error computation
        tic
        ur = u(:,1:i);
        [pe,~] = subsetselection(ur, 'pqr');
        Fh = ur*(ur(pe,:)\F(pe,:));
        svdtime_ = svdtime_ + toc;
        
        Dnrm(i) = norm(pinv(ur(pe,:)));
        
        Err = F-Fh;
        Errnrm(i) = Errnrm(i) + max(sqrt(sum(Err.^2,1)))./Fnrm;
        
        
        
        % RDEIM Error computation
        tic
        uh = rangefinder(S,i,10,1);
        [pr,~] = subsetselection(uh, 'pqr');
        Fr = uh*(uh(pr,:)\F(pr,:));
        rtime_ = rtime_ + toc;
        
        
        
        Err = F-Fr;
        Errnrmr(i) = Errnrmr(i) +  max(sqrt(sum(Err.^2,1)))./Fnrm; 
        
        Dhnrm(i) = Dhnrm(i) + norm(pinv(uh(pe,:)));
        count_intersect = count_intersect + length(intersect(pe,pr));
    end
    
    svdtime(i) = svdtime(i) + svdtime_/nreps;
    rtime(i) = rtime_/nreps;
    ratio_intersect(i) = count_intersect/(nreps*i);
    Dhnrm(i) = Dhnrm(i)/nreps;
end
Errnrm = Errnrm./nreps;
Errnrmr = Errnrmr./nreps;

%%
figure, semilogy(1:ns,Errnrm,':',1:ns,Errnrmr,'--','Linewidth',4)
legend('DEIM','RDEIM')
set(gca,'FontSize', 16)
title('Relative Error','FontSize',20)
xlabel('Basis size r','FontSize',18)
print 'figs/example1' -depsc

%% 
figure, semilogy(1:ns, svdtime, 'k-',1:ns,rtime, 'r--', 'LineWidth',2.0)
ylim([0.01,1.0])
legend({'DEIM','RDEIM'}, 'FontSize', 16)
set(gca,'FontSize', 16)
title('Computational time','FontSize',20)
xlabel('Basis size r','FontSize',18)
print 'figs/example1_timing' -depsc



