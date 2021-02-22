clc
close all


%%  Example 1a Rapidly decaying eigenvalues 2 norm
A = gallery('randsvd',1024,1.e10,3);
[q, errest, iter] = adaptrange(A,20,15,1.e-2);

b = q'*A;
disp('Example 1a')
sprintf('Error estimate %g, Actual Error %g, Iterations %d', errest, ...
                    norm(A - q*b)/norm(A), iter)
k = size(b,1);
figure, semilogy(1:size(A,1),svd(A),'k-', ...
            1:k,svd(b),'r--')
title('Example 1a')

%%  Example 1b Rapidly decaying eigenvalues F-norm
A = gallery('randsvd',1024,1.e10,3);
[q, errest, iter] = adaptrange(A,20,15,1.e-2, 'fro');

b = q'*A;
disp('Example 1b')
sprintf('Error estimate %g, Actual Error %g, Iterations %d', errest, ...
                    norm(A - q*b, 'fro')/norm(A,'fro'), iter)
k = size(b,1);
figure, semilogy(1:size(A,1),svd(A),'k-', ...
            1:k,svd(b),'r--')
title('Example 1b')        



%%  Example 1c Rapidly decaying eigenvalues F-norm
A = gallery('randsvd',1024,1.e10,3);
[q, errest, iter] = adaptrange_frob(A,20,15,1.e-2);

b = q'*A;
disp('Example 1c')
sprintf('Error estimate %g, Actual Error %g, Iterations %d', errest, ...
                    norm(A - q*b, 'fro')/norm(A,'fro'), iter)
k = size(b,1);
figure, semilogy(1:size(A,1),svd(A),'k-', ...
            1:k,svd(b),'r--')
title('Example 1c')        
        

        
%% Example 2 - Singular value gap 
% (modified from Sorensen and Embree, SISC, 2016)
drop = 1.e-3;
A = sparse(3000, 300); 
for k = 1:300
    x = sprand(3000, 1, 0.025); y = sprand(300, 1, 0.025);
    if k < 10
        A = A + (1/k)*x*y';
    else
        A = A + (drop/k)*x*y';
    end
    
end

[q, errest, iter] = adaptrange_frob(A,20,15,1.e-2);
b = q'*A;
disp('Example 2')
sprintf('Error estimate %g, Actual Error %g', errest, norm(A - q*b))
k = size(b,1);
figure, semilogy(1:size(A,2),svd(full(A)),'k-', ...
            1:k,svd(b),'r--')
title('Example 2')
        