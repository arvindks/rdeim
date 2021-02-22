clc
close all

%% Simulate the full order system
N = 1000;
x0 = zeros(N,1);

tspan = linspace(0,7,2000);
[tf,xf] = ode45(@(t,x) nonlinearrc(t,x, N), tspan, x0);




%%  Compute POD basis
[u,s,~] = svd(xf',0);
s = diag(s);
ind = find(cumsum(s.^2)/sum(s.^2) < 1-1.e-6);
V = u(:,ind);


%% Compute DEIM basis
nonlin = @(x) exp(40.0*x)+x-1;
A1 = spdiags(ones(N-1,1),-1,N,N)-speye(N);
A2 = spdiags([ones(N-1,1);0],0,N,N)-spdiags(ones(N,1),1,N,N);
fxf = nonlin(A1*xf')-nonlin(A2*xf');    
[W, ~, ~] = randQB_EI_auto(fxf, 1.e-3, 5, 1);

%% Compute optimal interpolating points using pivoted qr
%[~,~,p] = qr(W',0);
%k = size(ind,1);
%p = p(1:k);
[p,~,~] = hybrid(W);

%% Simulate reduced order model
VtW = V'*W;
x0 = zeros(size(V,2),1);
[tr,xh] = ode45(@(t,x) nonlinearrcrom(t,x, N, p, V, W, VtW), tspan, x0);

%% Compare ROM vs FOM
figure;
xr = V*xh';
subplot(1,2,1)
plot(tf,xf(:,1),'k-',tr,xh(:,1),'r--', 'LineWidth', 3.0)
legend('Full','ROM')
subplot(1,2,2)
plot(tf,abs(xf(:,1)-xh(:,1)),'r--', 'LineWidth', 3.0)
set(gca, 'FontSize', 18)