function [st] = angle_bounds(V,Omega,s,k,q)

    Omega  = V'*Omega;
    Omega1 = Omega(1:k,:);  Omega2 = Omega(k+1:end,:);
    Onrm = norm(Omega2/Omega1);
    
    s1 = s(k);    s2 = s(k+1);    
    
    tau = (s2./s1).^(2*q+1); 
    tO  = tau*Onrm;
    
    st = tO./sqrt(1+tO.^2);
        
end