function [N,Nx,Nxx,parm] = eqNcycle(parm,in,x)
% in is the mapping from x to parameter names (see switch below)
% N: nitrogen model state
% Fx: partial derivative of nitrogen model eqn w.r.t. model parameters,x
% Fn: partial derivative of nitrogen model eqn w.r.t. model state,N
% Fp: partial derivative of nitrogen model eqn w.r.t. p-model state,P
% Fxx: hessian matrix of N model wir.t. model parameters x
    
    global GN 
    iwet = find(parm.M3d(:));
    n_wet = length(iwet);
    % unpack the parameters to be optimized
    for ik1 = 1:length(in)
        switch in(ik1)
          case 1
            parm.bp = exp(x(ik1));
          case 2
            parm.kappa_dp = exp(x(ik1));
          case 3
            parm.alpha = exp(x(ik1));
          case 4
            parm.beta = exp(x(ik1));
          case 5
            parm.kw = exp(x(ik1));
          case 6
            parm.bn = exp(x(ik1));
          case 7
            parm.kappa_dn = exp(x(ik1));
          case 8
            parm.cc = exp(x(ik1));
          case 9
            parm.n2p_l = exp(x(ik1));
          case 10
            parm.S = exp(x(ik1));
          case 11
            parm.Nc = x(ik1);
        end
    end
    
    % use the global variable to set the initial iterate for the
    % N-cycle Newton solver
    DIN = GN(1:n_wet);
    PON = GN(n_wet+1:2*n_wet);
    DON = GN(2*n_wet+1:end);
    X0  = [DIN;PON;DON];
    
    %
    % Solve the N-cycle equilibrium 
    %
    options.atol = 5e-15; options.rtol = 1e-16; options.iprint = 1;
    [N,ierr]...
        = nsnew(X0,@(X) N_eqn(X,parm,in,x),options);
    if (ierr ~=0)
        fprintf('eqNcycle did not converge.\n');
        keyboard
    else
        GN = 0.999*real(N); % reset the global variable for the next call eqNcycle
        
        if nargout>1     
            %
            % Compute the gradient of the solution wrt the parameters
            %
            [F,FFn,Fn,Nx,Nxx,parm] = N_eqn(N,parm,in,x);
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,FFn,Fn,Nx,Nxx,parm] = N_eqn(X,parm,in,x)
% unpack some useful stuff
    c2n = parm.c2n;
    
    M3d = parm.M3d;
    TRdiv = parm.TRdiv;
    if (~isempty(find(in==1)))
        % do this if b needs to be optimized
        [PFdiv_n,dPFDdb,d2PFDdb2,dPFDdkappa_p] ...
            = buildPFD(M3d,parm.bn,parm.kappa_p,parm.grd);
    else
        % do this if the  parameters don't need to be optimized
        PFdiv_n = buildPFD(M3d,parm.bn,parm.kappa_p,parm.grd);  % particle flux divergence [s^-1];
    end
    DIP = parm.DIP;
    
    L = parm.L;
    Pass = parm.alpha*(L*DIP);
    
    iwet   = find(M3d(:));
    n_wet = length(iwet);
    I      = speye(n_wet);
    
    DIN           = X(1:n_wet);
    PON           = X(n_wet+1:2*n_wet);
    DON           = X(2*n_wet+1:end);
    POC           = PON.*parm.c2n;
    % make an N:P function
    % parm.n2p   = M3d(iwet)*0;
    % mm = find(parm.po4obs(iwet) >= parm.DIPc);
    % parm.n2p(mm) = parm.n2p_l;
    % nn = find(parm.po4obs(iwet) < parm.DIPc);
    % parm.n2p(nn) = parm.n2p_l-parm.S*(parm.po4obs(iwet(nn))-parm.DIPc);
    parm.n2p = parm.n2p_l+parm.S*(1-tanh(parm.po4obs(iwet)-parm.DIPc));
    
    % make a water column denitrification function
    parm.io = 0.5-0.5*tanh((parm.o2(iwet)-parm.Oc)/parm.cc);
    %IX      =  parm.sl*0.99.^(parm.o2(iwet)-DIN)+parm.iN;
    %dIXdDIN = -parm.sl*0.99.^(parm.o2(iwet)-DIN)*log(0.99);
    %d2IXdDIN2 = -parm.sl*0.99;
    %SDN     = -(PFdiv_n*POC).*parm.ib.*IX;
    
    
    [G,Gx,Gxx] = uptake(parm);
    [wcD,sdD,Dn,rhsDx] = denit_don(X,PFdiv_n,dPFDdb,d2PFDdb2,parm);
    [R,Rn,rhsRx] = fixit(X,parm);
    %WARNING
    % D = k_DN.*POC.*DIN+SDN
    F = ([TRdiv*DIN+R.*G-parm.kappa_dn*DON+wcD+sdD-parm.DIN_dep-parm.DIN_RIV;...
          -(1-parm.sigma)*G+(PFdiv_n+parm.kappa_p*I)*PON;...
          -parm.sigma*G+(TRdiv+parm.kappa_dn*I)*DON-parm.kappa_p*PON-parm.DON_RIV]);
    
    if nargout>1
        % Jacobian w.r.t. n-model state
        Fn = [[TRdiv+d0(G.*Rn),          0*I,-parm.kappa_dn*I]+Dn;...
              [ 0*I,  PFdiv_n+parm.kappa_p*I,             0*I];...
              [ 0*I,         -parm.kappa_p*I, TRdiv+parm.kappa_dn*I]];
        
        % Jacobian w.r.t. DIP
        Fp = [[          d0(R)*parm.alpha*L*(d0(parm.n2p)), 0*I,0*I];...
              [-(1-parm.sigma)*parm.alpha*L*(d0(parm.n2p)), 0*I,0*I];...
              [-(  parm.sigma)*parm.alpha*L*(d0(parm.n2p)), 0*I,0*I]];
        
        % Factored Jacobian
        FFn = mfactor(Fn);
        
        
    end    
    if nargout>3
        %DIN = N(1:n_wet);
        %PON = N(n_wet+1:2*n_wet);
        %DON = N(2*n_wet+1:end);
        
        DIPx = parm.Px(1:n_wet,:);
        Z = zeros(n_wet,1);
        %  WARNING
        % for case 1-4 use real to test Nx not for testing Nxx
        for ik1 = 1:length(in)
            switch in(ik1)
              case 1 %  bp
                Fx(:,ik1) = [R.*Gx(:,1);...
                             -(1-parm.sigma)*Gx(:,1);...
                             -(parm.sigma)*Gx(:,1)];
              case 2 % kappa_dp
                Fx(:,ik1) = [R.*Gx(:,2);...
                             -(1-parm.sigma)*Gx(:,2);...
                             -(parm.sigma)*Gx(:,2)];
              case 3  % alpha
                      %dPassdalpha = parm.L*parm.DIP;
                      %Fx(:,ik1) = exp(x(ik1))*[R.*dPassdalpha.*parm.n2p;...
                      %       -(1-parm.sigma)*dPassdalpha.*parm.n2p;...
                      %       -parm.sigma*dPassdalpha.*parm.n2p]+Fp*parm.Px(:,3);
                
                Fx(:,ik1) = [R.*Gx(:,3);...
                             -(1-parm.sigma)*Gx(:,3);...
                             -(parm.sigma)*Gx(:,3)];
              case 4  % beta
                      %dPassdbeta = parm.alpha*parm.dLdbeta*DIP;
                      %Fx(:,ik1) = exp(x(ik1))*[R.*dPassdbeta.*parm.n2p;...
                      %       -(1-parm.sigma)*dPassdbeta.*parm.n2p;...
                      %       -parm.sigma*dPassdbeta.*parm.n2p]+Fp*parm.Px(:,4);
                Fx(:,ik1) = [R.*Gx(:,4);...
                             -(1-parm.sigma)*Gx(:,4);...
                             -(parm.sigma)*Gx(:,4)];
                
              case 5  % kw 
                
                %Fx(:,ik1) = exp(x(ik1))*[parm.io.*POC.*DIN;...
                %             Z;...
                %             Z];
                Fx(:,ik1) = [rhsDx(:,ik1);Z;Z];
              case 6  % bn
                
                %Fx(:,ik1) = exp(x(ik1))*[-(dPFDdb*POC).*parm.ib.*IX;...
                %              (dPFDdb*PON);...
                %              Z];
                Fx(:,ik1) = [rhsDx(:,ik1);...
                             exp(x(ik1))*(dPFDdb*PON);...
                             Z];
                %              case 2  % sigma
                % Fx(:,ik1) = [Z;...
                %             -Pass.*parm.n2p;
                %             -Pass.*parm.n2p];
              case 7  % kappa_dn
                Fx(:,ik1) = exp(x(ik1))*[-DON;...
                                    Z;...
                                    DON];
              % case 8  % Oc
              
                %diodOc = -(tanh((parm.Oc-parm.o2(iwet))/parm.cc).^2-1)/ ...
                %         (2*parm.cc);
                %Fx(:,ik1) = exp(x(ik1))*[parm.kw*diodOc.*POC.*DIN;...
                %                    Z;...
                %                    Z];
                % Fx(:,ik1) = [rhsDx(:,ik1);Z;Z];
              case 8  % cc
                
                %diodcc = ((parm.Oc-parm.o2(iwet)).*...
                      %    (tanh((parm.Oc-parm.o2(iwet))/parm.cc).^2 - 1))/(2*parm.cc^2);
                
                %Fx(:,ik1) = exp(x(ik1))*[parm.kw*diodcc.*POC.*DIN;...
                %                    Z;...
                %                    Z];
                Fx(:,ik1) = [rhsDx(:,ik1);Z;Z];
              case 9  % n2p_l
                       %dn2pdn2p_l = 1;
                       %parm.dn2pdn2p_l = 1;
                       %Fx(:,ik1) = exp(x(ik1))*[ R.*Pass.*dn2pdn2p_l;...
                       %       -(1-parm.sigma)*Pass.*parm.dn2pdn2p_l;...
                       %       -parm.sigma*Pass.*parm.dn2pdn2p_l];
                Fx(:,ik1) = [R.*Gx(:,9);...
                             -(1-parm.sigma)*Gx(:,9);...
                             -parm.sigma*Gx(:,9)];
              case 10 % S
                      %dn2pdS     = M3d(iwet)*0;
                      %dn2pdS(mm)     = 0;
                      %dn2pdS(nn)     = -(parm.po4obs(iwet(nn))-parm.DIPc);    
                      %parm.dn2pdS = dn2pdS;
                      %Fx(:,ik1) = exp(x(ik1))*[ R.*Pass.*parm.dn2pdS;...
                      %        -(1-parm.sigma)*Pass.*parm.dn2pdS;...
                      %        -parm.sigma*Pass.*parm.dn2pdS];
                Fx(:,ik1) = [R.*Gx(:,10);...
                             -(1-parm.sigma)*Gx(:,10);...
                             -parm.sigma*Gx(:,10)];
              case 11 % Nc
                %[~,~,rhsRx] = fixit(X,parm);
                %dTdNc  = -(1-T.^2);
                %dRdNc  = 0.5*dTdNc;
                Fx(:,ik1) = [G.*rhsRx;...
                             Z;...
                             Z];
            end
        end
        Nx = mfactor(FFn,-Fx);
        
    end
    if (nargout>5)
        dN = mfactor(FFn,-F);
        N = X+dN;
        
        nin = parm.nin;
        pos = parm.pos;
        ZZ = zeros(n_wet,pos(end,end));
        rhsFxx = zeros(3*n_wet,pos(end,end));
        %
        DINx = Nx(1:n_wet,:);
        PONx = Nx(n_wet+1:2*n_wet,:);
        DONx = Nx(2*n_wet+1:end,:);
        %WARNING
        [~,~,~,Rx,rhsRxx] = fixit(N,parm,Nx);
        % WARNING
        [~,~,~,~,Dx,rhsDxx] = denit_don(N,PFdiv_n,dPFDdb,d2PFDdb2,parm,Nx);
        %$ Rx Gx 
        for i1 = 1:nin
            for i2 = i1:nin
                rhsFxx(:,pos(i1,i2)) = rhsFxx(:,pos(i1,i2))+...
                    [Rx(:,i1).*Gx(:,i2)+Rx(:,i2).*Gx(:,i1)+...
                     R.*Gxx(:,pos(i1,i2))+ G.*rhsRxx(:,pos(i1,i2))+rhsDxx(:,pos(i1,i2));...
                     -(1-parm.sigma)*Gxx(:,pos(i1,i2));...
                     -parm.sigma*Gxx(:,pos(i1,i2))];
            end
        end
        % bn do an L
        % vertical part of L
        for i1 = 1:6
            rhsFxx(:,pos(i1,6)) = rhsFxx(:,pos(i1,6))+...
                   [ZZ(:,1);...
                    parm.bn*dPFDdb*PONx(:,i1);...
                    ZZ(:,1)];
        end
        % corner part of L
        rhsFxx(:,pos(6,6)) = rhsFxx(:,pos(6,6))+...
            [ZZ(:,1);...
             (parm.bn*dPFDdb+parm.bn^2*d2PFDdb2)*PON;...
             ZZ(:,1)];
        % horizontal part of L
        r = [pos(6,6):pos(6,end)];
        rhsFxx(:,r) = rhsFxx(:,r) +...
            [ZZ(:,r);...
             parm.bn*dPFDdb*PONx(:,6:end);...
             ZZ(:,r)];
        
        % kappa_dn
        gk = zeros(n_wet,11);
        gk(:,7) = ones(n_wet,1)*parm.kappa_dn;
        ggk = zeros(n_wet,pos(end,end));
        ggk(:,pos(7,7)) = ones(n_wet,1)*parm.kappa_dn;
        for i1= 1:11
            for i2 = i1:11
                tmp = (ggk(:,pos(i1,i2)).*DON + gk(:,i1).*DONx(:,i2)+gk(:,i2).*DONx(:,i1));
                rhsFxx(:,pos(i1,i2)) = rhsFxx(:,pos(i1,i2))+...
                    [-tmp;...
                     Z;...
                     tmp];
            end
        end
        
        %
        Nxx = mfactor(FFn,-rhsFxx);
    end
    
