function [P,Px,Pxx,parm] = eqPcycle(parm,ip,x)
% ip is the mapping from x to parameter names (see switch below)
% output: P is model prediction of DIP,POP,and DOP
% output: F partial derivative of P model w.r.t. model parameters x
% output: Fxx hessian matrix of P model w.r.t.  model parameters x
    
% unpack the parameters to be optimized
    if (nargin>1)
        for ik1 = 1:length(ip)
            switch ip(ik1)
              case 1
                parm.bp = exp(x(ik1)); % Martin exponent for POP solubilization
              case 2
                parm.sigma = exp(x(ik1)); % fraction of organic P allocated to dissolved pool
              case 3
                parm.kappa_dp = exp(x(ik1)); % DOP remineralization const.[s^-1];
              case 4
                parm.alpha = exp(x(ik1)); % npp scaling factor for DIP uptake rate
              case 5
                parm.beta = exp(x(ik1)); % npp scaling exponent for DIP uptake rate
            end
        end
    end
    bp = parm.bp;             
    sigma = parm.sigma;       
    kappa_dp = parm.kappa_dp; 
    alpha = parm.alpha;       
    beta = parm.beta;         
    
    % unpack some useful stuff
    M3d = parm.M3d;
    TRdiv = parm.TRdiv;
    grd = parm.grd;
    
    iwet = find(M3d(:));         % wet point indices;
    n_wet = length(iwet);        % number of wet points;
    
    I = speye(n_wet);            % make an identity matrix;
    d0 = @(x) spdiags([x(:)],[0],length(x(:)),length(x(:)));
    
    % fixed parameters
    DIPbar = M3d(iwet)*parm.DIPbar;  % gobal arerage PO4 conc.[mmol m^-3]; 
    kappa_g = parm.kappa_g;          % PO4 geological restore const.[s^-1];
    kappa_p = parm.kappa_p;          % POP solubilization rate constant
    npp = parm.npp;                  % net primary production 
    
    % build part of the biological DIP uptake operator
    Lambda = parm.Lambda;
    LAM = 0*M3d;
    LAM(:,:,1) = (npp.^beta).*Lambda(:,:,1);
    LAM(:,:,2) = (npp.^beta).*Lambda(:,:,2);
    L = d0(LAM(iwet));  % PO4 assimilation rate [s^-1];
    parm.L = L; 
    
    
    % build the sinking particulate flux divergence operator
    if (~isempty(find(ip==1)))
        % do this if b needs to be optimized
        [PFdiv,dPFDdb,d2PFDdb2] ...
            = buildPFD(M3d,bp,kappa_p,grd);
    else
        % do this if the  parameters don't need to be optimized
        PFdiv = buildPFD(M3d,bp,kappa_p,grd);  % particle flux divergence [s^-1];
    end
    
    % build Jacobian matrix
    %    [                 dF/dDIP,          dF/dPOP,          dF/dDOP]            
    Fp = [ TRdiv+alpha*L+kappa_g*I,              0*I,       -kappa_dp*I; ...  % DIP eqn
                -(1-sigma)*alpha*L,  PFdiv+kappa_p*I,               0*I; ...  % POP eqn            
                    -sigma*alpha*L,       -kappa_p*I,  TRdiv+kappa_dp*I];     % DOP eqn
    
    % right hand side of phosphate equations 
    Z=zeros(n_wet,1);
    RHS = [ kappa_g*DIPbar;...
                       Z;...
                       Z];
    
    % dP/dt + Fp*P = RHS    
    % factoring Jacobian matrix
    FFp = mfactor(Fp); 
    % solve for P-cycle model state
    P = mfactor(FFp,RHS);    
    if (nargout>1)
        %
        %
        % Compute the gradient of the solution wrt the parameters
        %
        %
        DIP = P(1:n_wet);
        POP = P(n_wet+1:2*n_wet);
        DOP = P(2*n_wet+1:end);
        Fx = zeros(3*n_wet,length(ip));
        for ik1 = 1:length(ip)
            switch ip(ik1)
              case 1 % b                
                Fx(:,ik1) = exp(x(ik1))* ...
                    [Z;...
                     dPFDdb*POP;...
                     Z];
              case 2 % sigma
                Fx(:,ik1) =  exp(x(ik1))* ...
                    [Z;...
                     alpha*L*DIP;...
                     -alpha*L*DIP];
              case 3 % kappa_dp
                Fx(:,ik1) = exp(x(ik1))* ...
                    [-DOP;...
                     Z;...
                     DOP];
              case 4 % alpha
                Fx(:,ik1) = exp(x(ik1))* ...
                    [L*DIP;...
                     -(1-sigma)*L*DIP;...
                     -sigma*L*DIP];
                
                
              case 5 %beta
                dLambdadbeta = 0*Lambda;
                dLambdadbeta(:,:,1) = log(npp).*LAM(:,:,1);
                dLambdadbeta(:,:,2) = log(npp).*LAM(:,:,2);
                iz = find(isinf(dLambdadbeta(:)));
                dLambdadbeta(iz) = 0;
                inan = find(isnan(dLambdadbeta(:)));
                dLambdadbeta(inan) = 0;
                dLdbeta = d0(dLambdadbeta(iwet)); 
                Fx(:,ik1) = exp(x(ik1))*[ alpha*dLdbeta*DIP;...
                                    -(1-sigma)*alpha*dLdbeta*DIP;...
                                    -sigma*alpha*dLdbeta*DIP];
                % will need L and dLdbeta for gradients of other
                % biogeochemical cycles
                parm.L = L; 
                parm.dLdbeta = dLdbeta;
            end
        end
        % Fx is the derivative of the solution wrt to the parameters
        Px = mfactor(FFp,-Fx);
    end
    if (nargout > 2)
        %
        %
        % Compute the 2nd derivative of the solution wrt the parameters
        %
        %
        
        % initialize the 2nd derivative of the Jacobian wrt the parameters
        for rr = 1:3 % 3 is the number of species DIP,POP,DOP
            for cc = 1:3 % 3 is the number of species DIP,POP,DOP
                for i1 = 1:5 % 5 is the number of parameters, i.e. b,kappa_dp,sigma,alpha,beta
                    for i2 = i1:5 % 5 is the number of parameters
                        d2J{rr,cc,i1,i2} = sparse(n_wet,n_wet);
                    end
                end
            end
        end
        DIPx = Px(1:n_wet,:);
        POPx = Px(n_wet+1:2*n_wet,:);
        DOPx = Px(2*n_wet+1:end,:);
        % compute only the upper triangular part of the matrix
        Z = zeros(n_wet,length(ip));
        for i1 = 1:length(ip)
            switch ip(i1)
              case 1 % b
                FpxPx(:,:,i1) = exp(x(i1))*[  Z;...
                                  dPFDdb*POPx;...
                                  Z];
                for i2 = i1:length(ip) 
                    switch ip(i2)
                      case 1 % b,b
                        d2J{2,2,1,1} = d2PFDdb2;
                      case 2 % b,sigma
                      case 3 % b,kappa_dp
                      case 4 % b,alpha
                      case 5 % b,beta
                    end
                end
              case 2 % sigma
                FpxPx(:,:,i1) =  exp(x(i1))*[  Z;...
                                    alpha*L*DIPx;...
                                   -alpha*L*DIPx];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 2 % sigma,sigma
                      case 3 % sigma,kappa_dp
                      case 4 % sigma,alpha
                        d2J{2,1,2,4} =  L;
                        d2J{3,1,2,4} = -L;
                      case 5 % sigma,beta
                        d2J{2,1,2,5} =  alpha*dLdbeta;
                        d2J{3,1,2,5} = -alpha*dLdbeta;
                    end
                end
              case 3 % kappa_dp
                FpxPx(:,:,i1) = exp(x(i1))*[  -DOPx;...
                                  Z;...
                                  DOPx];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 3 % kappa_dp,kappa_dp
                      case 4 % kappa_dp,alpha
                      case 5 % kappa_dp,beta
                    end
                end
              case 4 % alpha
                FpxPx(:,:,i1) = exp(x(i1))*[  L*DIPx;...
                                 -(1-sigma)*L*DIPx;...
                                 -sigma*L*DIPx];
                
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 4 % alpha,alpha
                      case 5 % alpha,beta
                        d2J{1,1,4,5} = dLdbeta;
                        d2J{2,1,4,5} = -(1-sigma)*dLdbeta;
                        d2J{3,1,4,5} =     -sigma*dLdbeta;
                    end
                end
              case 5 % beta
                FpxPx(:,:,i1) = exp(x(i1))*[  alpha*dLdbeta*DIPx;...
                                  -(1-sigma)*alpha*dLdbeta*DIPx;...
                                  -sigma*alpha*dLdbeta*DIPx];
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 5 % beta,beta
                        d2Lambdadbetadbeta = 0*Lambda;
                        d2Lambdadbetadbeta(:,:,1) = log(npp).*log(npp).*LAM(:,:,1);
                        d2Lambdadbetadbeta(:,:,2) = log(npp).*log(npp).*LAM(:,:,2);
                        iz = find(isinf(d2Lambdadbetadbeta(:)));
                        d2Lambdadbetadbeta(iz) = 0;
                        inan = find(isnan(d2Lambdadbetadbeta(:)));
                        d2Lambdadbetadbeta(inan) = 0;
                        d2Ldbetadbeta = d0(d2Lambdadbetadbeta(iwet)); 
                        d2J{1,1,5,5} = alpha*d2Ldbetadbeta;
                        d2J{2,1,5,5} = -(1-sigma)*alpha*d2Ldbetadbeta;
                        d2J{3,1,5,5} =     -sigma*alpha*d2Ldbetadbeta;
                        parm.d2Ldbetadbeta = d2Ldbetadbeta;
                    end
                end
            end
        end
        symcol = @(i,j) (i>=j).*i+(i<j).*j;
        symrow = @(i,j) (i>=j).*j+(i<j).*i;
        delta = @(i,j) (i==j)*1 + (i~=j)*0.0;
        k = 0;
        for i1 = 1:length(ip)
            for i2 = i1:length(ip) 
                k = k+1;
                rhs(:,k) = -(...
                    exp(x(i1)+x(i2))*cell2mat(d2J(:,:,symrow(ip(i1),ip(i2)),symcol(ip(i1),ip(i2))))*P+...
                    FpxPx(:,i1,i2)+FpxPx(:,i2,i1)+delta(i1,i2)*Fx(:,i1));

            end
        end
        Pxx = mfactor(FFp,rhs);
    end
    