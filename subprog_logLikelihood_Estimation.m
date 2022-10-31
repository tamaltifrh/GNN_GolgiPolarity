function logL=subprog_logLikelihood_Estimation(mu,A,B,g)    

    n = size(A,1); m = size(A,2); 
    b = zeros(n,1);
    tau = sqrt(mu);
    
    S1 = [A      b
          tau*B  tau*g];
  
      if exist('spqr')==3 && size(S1,2)>2000  %if SuiteSparse SPQR is installed
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
    %Stild=Qtild*Rtild; %Rtild, Qtild : new QR decomposition with R diagonal positive
    
    H = Rtild(1:m,1:m);
    h_bold = Rtild(1:m,m+1);
    h=Rtild(m+1,m+1);
    
    E_mu = h^2;
    
    BB=B'*B; % B'B is diagonal with 1s and 0s
    diagBB=diag(BB);
    m0=length(diagBB)-sum(diagBB); % number of zero eigenvalues of B'B
    
    M=mu*BB;
    diagM=diag(M);
    diagM(diagM==0)=[];
    
    diagH=diag(H(1:m-1,1:m-1));
    if length(find(diagH))~=length(diagH)
        diagH(diagH==0)=[]; % beware of null diagonal elements for the log (there shouldn't be any)
        disp('Beware, segmentation problem: Null diagonal elements of H had to be removed. Hunt down white cells in the TissueAnalyzerToMatlab outcome !');
    end   
    
    logL=-(n-m0+1)*log(E_mu)+sum(log(diagM))-2*sum(log(diagH));
 
end