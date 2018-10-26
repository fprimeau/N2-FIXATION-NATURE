function [G,Gx,Gxx] = uptake(parm)
% unpack the parameters to be optimized
alpha = parm.alpha;
beta = parm.beta;
S = parm.S;
n2p_l = parm.n2p_l;
dLdbeta = parm.dLdbeta;
d2Ldbetadbeta = parm.d2Ldbetadbeta;
M3d = parm.M3d; iwet = find(M3d(:)); n_wet = length(iwet);
% N:P of uptake operator
po4obs = parm.po4obs(iwet);
n2p = parm.n2p_l+parm.S*(1-tanh(po4obs-parm.DIPc));
% DIPc = parm.DIPc;
%n2p = M3d(iwet)*0;
% mm = find(po4obs(iwet)>=parm.DIPc); 
% nn = find(po4obs(iwet)<parm.DIPc);
% n2p(mm) = n2p_l;
% n2p(nn) = n2p_l-S*(po4obs(iwet(nn))-DIPc);
dn2pdn2p_l = 1;
dn2pdS = 1-tanh(po4obs-parm.DIPc);
% dn2pdS(nn) = -(po4obs(iwet(nn))-DIPc);
% N uptake operator
L = diag(parm.L);  
dLdbeta = diag(parm.dLdbeta); 
d2Ldbetadbeta = diag(parm.d2Ldbetadbeta);

DIP = parm.DIP;
G = alpha.*n2p.*L.*DIP; 
Gp = alpha*n2p.*L;   

% gradient of uptake operator
nin = 11;
Gpx = zeros(n_wet,nin);
Gpx(:,3) = alpha.*n2p.*L;              % dGdlog_alpha
Gpx(:,4) = beta*alpha.*n2p.*dLdbeta;   % dGdlog_beta
Gpx(:,9) = n2p_l*alpha.*dn2pdn2p_l.*L; % dGdlog_n2p_l
Gpx(:,10) = S*alpha*dn2pdS.*L;          % dGdlog_S

% Gradient
% grad DIP
DIPx = zeros(n_wet,11);
DIPx(:,1:4) = parm.Px(1:n_wet,:);   

Gx = zeros(n_wet,11);
Gx = d0(Gp)*DIPx+d0(DIP)*Gpx;

% Hessian
% map P-model parameter order to N-model parameter order
nin = parm.nin;
n = (nin+1)*nin/2;
pos = tril(ones(nin),0);
pos(pos==1) = 1:n;
pos = pos';

% grad grad DIP
DIPxx = zeros(n_wet,n);
map = [pos(1,1) pos(1,2) pos(1,3) pos(1,4) pos(2,2) pos(2,3) pos(2,4) ...
       pos(3,3) pos(3,4) pos(4,4)];
DIPxx(:,map) = parm.Pxx(1:n_wet,:);

Gpxx = zeros(n_wet,n);
Gpxx(:,pos(3,3)) = alpha*n2p.*L;                          % alpha alpha
Gpxx(:,pos(3,4)) = beta*alpha*n2p.*dLdbeta;               % alpha beta
Gpxx(:,pos(3,9)) = n2p_l*alpha*dn2pdn2p_l.*L;            % alpha n2p_l 
Gpxx(:,pos(3,10)) = S*alpha*dn2pdS.*L;                    % alpha S 
Gpxx(:,pos(4,4)) = beta*alpha*n2p.*dLdbeta+...          
    beta^2*alpha*n2p.*d2Ldbetadbeta;                      % beta beta 
Gpxx(:,pos(4,9)) = n2p_l*beta*alpha*dn2pdn2p_l.*dLdbeta; % beta n2p_l
Gpxx(:,pos(4,10)) = S*beta*alpha*dn2pdS.*dLdbeta;         % beta S
Gpxx(:,pos(9,9)) = n2p_l*alpha*dn2pdn2p_l.*L;           % n2p_l n2p_l
Gpxx(:,pos(9,10)) = 0*G;                                 % n2p_l S    
Gpxx(:,pos(10,10)) = S*alpha*dn2pdS.*L;                   % S S

%
Gxx = d0(DIP)*Gpxx+d0(Gp)*DIPxx;
for is = 1:nin
    r = [pos(is,is):pos(is,end)];
    Gxx(:,r) = Gxx(:,r)+...
        d0(Gpx(:,is))*DIPx(:,is:end)+...
        d0(DIPx(:,is))*Gpx(:,is:end);
end





