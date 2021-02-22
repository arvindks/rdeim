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



%% 


%% Test angle bounds
[u,s,v] = svd(S,0);
s = diag(s);
n = nx*ny;
theta = zeros(4,5);
thetab = zeros(4,5);

p = 20;
stbnd = @(r,p,n,g,q) (sqrt(r/(p-1))+exp(1)*sqrt((r+p)*(n-r))/p)*g^(2*q+1)/(1-g); 

   i = 0;
for k = 30:20:90
    i = i +1 ;  g = s(k+1)/s(k);
    Omega = randn(n, k + p);    
    
    vo = v'*Omega;
    O1 = vo(1:k,:); O2 = vo(k+1:end,:);
    
    for q = 0:4
        Q = subspaceiteration(S, Omega, q);
        [uh,~,~] = svd(Q'*S,'econ');
        
        tb = subspace_angles(u(:,1:k),Q*uh(:,1:k));
        theta(i,q+1) = sin(max(tb));
        
        bnd = norm(O2*pinv(O1))*g.^(2*q+1)/(1-g);
        %bnd = stbnd(k,p,n,g,q);
        thetab(i,q+1) = bnd;
        
    end
end


figure, semilogy(theta')
figure, semilogy(thetab')


%%
thetaq = zeros(4,20); i = 0;
for k = 30:20:90
    i = i +1 ;  g = s(k+1)/s(k);
    vo = v'*Omega;
    O1 = vo(1:k,:); O2 = vo(k+1:end,:);
    for q = 0:19
        
        %bnd = stbnd(k,p,n,g,q);
        Onrm = norm(O2/O1);
        bnd = Onrm*g^(2*q+1)/sqrt(1+ Onrm^2*g^(4*q+2)) ;
        thetaq(i,q+1) = bnd;
        
    end
end
 
figure, semilogy(thetaq')
