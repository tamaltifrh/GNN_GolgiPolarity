function [INFERENCE]=subprog_inverse_matrix_A(mu, A, B, g, E_tot)
    
    n = size(A,1); m = size(A,2); 
    b = zeros(n,1);
    tau = sqrt(mu);
    
    S1 = [A      b
          tau*B  tau*g];
    
      if exist('spqr')==3 %if SuiteSparse SPQR is installed (http://faculty.cse.tamu.edu/davis/suitesparse.html)
          S1=sparse(S1);
          R=spqr(S1);
      else
          [~,R]=qr(S1);
      end
    
    %matrix of signs of diagonal elements of R
    SR=zeros(m+n);
    for i=1:m+1
        SR(i,i)=sign(R(i,i));
    end
    %SRm1=SR; %this matrix is its own inverse
    Rtild=SR*R;
    %Qtild=Q*SRm1;
    %Stild=Qtild*Rtild; %Rtild, Qtild : new QR decomposition with R diagonal positive. 
    %Stild should be equal to S
    
    H = Rtild(1:m,1:m);
    H=sparse(H);
    Hm1=pseudoinverse(H);
    %Hm1=pinv(H);
    h_bold = Rtild(1:m,m+1);
    p = Hm1*h_bold;
    INFERENCE.TENSIONS  = p(1:E_tot)';
    INFERENCE.PRESSURES = p(E_tot+1:end)'-mean(p(E_tot+1:end)); % pressure is at zero average
    
end