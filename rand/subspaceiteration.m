function Q = subspaceiteration(A, Omega, q)

    Y = A*Omega;    [Q,~] = qr(Y,0);
    
    for i = 1:q
        y1 = A'*Q;   [Q,~] = qr(y1,0);
        y2 = A*Q;    [Q,~] = qr(y2,0);   
    end

end