clc 
close all

clc
close all

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

%Collect snapshots
S = zeros(nx*ny,size(mu,1));
for i = 1:size(mu,1)
   sv = s(X,Y,mu(i,:));
   S(:,i) =  sv(:);
end


n = nx*ny;

%% Test angle bounds
[u,s,v] = svd(S,0);

k = 100; p = 30; ell = k + p;
theta = zeros(k,1);
thetab = zeros(30,1); thetab2 = thetab1;

stbnd = @(r,p,skp1,sf) sqrt(skp1.^2+r*sf.^2/(p-1));
stbnd2 = @(r,p,n,g) (sqrt(r/(p-1))+exp(1)*sqrt((r+p)*(n-r))/p)*g; 
s = diag(s);
omega = randn(size(S,2),ell);
for i = 1:k
    Y = S*omega(:,1:i+p); [q,~] = qr(Y,0);   [uh,~,~] = svd(q'*S); ur = q*uh(:,1:i);
    theta(i) = sin(max(subspace_angles(u(:,1:i),ur)));
    %[st] = angle_bounds(v,omega,diag(s),i,0);
    

    g = s(i+1)/s(i);
    thetab(i) = stbnd2(i,p,n,g);
    
end

figure, semilogy(1:k,theta,1:k,thetab);
legend('real','bnd2')

