function [f,fx,fxx] = neglogpost(x,parm)
tic
M3d  = parm.M3d;
iwet = find(M3d(:));
nwet = length(iwet);

po4raw = parm.po4raw;
no3raw = parm.no3raw;
DONobs = parm.DONobs-1.8; % get rid of refractory DON
DONobs(:,:,6:end) = 0;    % use only upper 500 m DON for model constants.
scale = parm.DON_scale;
fprintf('DON scale is %2.2f, \n', scale)
fprintf('SDN scale is %2.2f, \n', parm.SDN_scale)

%%%%%%% check if the optimization routin gives strange parameter values
xii_new = [exp(x(1:end-1));x(end)];
load tmp_xhat_preind.mat
xii_old = xhat;
a = 5;

fprintf('current parameter scale is %2d, \n', a)

tmp = find(abs(xii_new)>=abs(a*xii_old) | ...
           abs(xii_new)<= abs((1/a)*xii_old));
for ii = 1:length(tmp)
    xii_new(tmp(ii)) = xii_old(tmp(ii));
end

% restrict b value below 1 and .3.
if (xii_new(1)>1.5 | xii_new(1)<0.5);
    xii_new(1) = xii_old(1);
end

if (xii_new(4)>1.25*xii_old(4) | ...
    xii_new(4)<0.75*xii_old(4));
    xii_new(4) = xii_old(4);
end

if (xii_new(6)>1.5 | xii_new(6)<0.5);
    xii_new(6) = xii_old(6);
end
for i1 = 1:length(xii_new)
    fprintf('%3.3e,',xii_new(i1));
end
fprintf('\n');

xhat = xii_new;
save('tmp_xhat_preind','xhat')
x = [log(xii_new(1:end-1));xii_new(end)];
%%%%%%% restore it back to old value %%%%%%%%%%%%%.

ipo4 = find(~isnan(po4raw(iwet))); % index for valid DIP measurements;
ino3 = find(~isnan(no3raw(iwet))); % index for valid DIN measurements;
idon = find(DONobs(iwet)>0);       % index for valid DON measurements;

n_po4 = length(ipo4); % number of valid PO4 obs.;
n_no3 = length(ino3); % number of valid NO3 obs.;
n_don = length(idon); % number of valid DON obs.;

% extract valide data;
NO3   = no3raw(iwet(ino3)); 
PO4   = po4raw(iwet(ipo4));
DON_o = DONobs(iwet(idon));
parm.DONstd = std(DON_o);

DIN = M3d+nan; 
PON = M3d+nan; 
DON = M3d+nan;

% get P model only parameters;
xp = x(1:4);
ip = [1,3,4,5];
[P,Px,Pxx,parm] = eqPcycle(parm,ip,xp); % P cycle model;
parm.P = P;     % vector for P field 
parm.Px = Px;   % first derivatives;
parm.Pxx = Pxx; % second derivatives;

DIP = M3d+nan;  DIP(iwet) = P(1:nwet);
POP = M3d+nan;  POP(iwet) = P(1+nwet:2*nwet);
DOP = M3d+nan;  DOP(iwet) = P(1+2*nwet:end);
parm.DIP = DIP(iwet); parm.POP = POP(iwet); parm.DOP = DOP(iwet);

% get N model parameters;
xn = x;
in = 1:length(x);
[N,Nx,Nxx,parm] = eqNcycle_v2(parm,in,xn);

DIN = M3d+nan;  DIN(iwet) = N(1:nwet);
PON = M3d+nan;  PON(iwet) = N(1+nwet:2*nwet);
DON = M3d+nan;  DON(iwet) = N(1+2*nwet:end);

% save N and P fields;
fname =sprintf('PN_DON%1.0e_SDN%d_preind',parm.DON_scale, parm.SDN_scale*100);
save(fname,'DIN','DON','PON','DIP','DOP','POP')
% covariance matrix is diagonal with variances given by
% the inverse of the grid box volume times the variance of the obs
dVt = parm.dVt;
Wp = d0(dVt(iwet(ipo4))./(parm.DIPstd.^2*sum(dVt(iwet))));
Wn = d0(dVt(iwet(ino3))./(parm.DINstd.^2*sum(dVt(iwet))));
Wo = d0(dVt(iwet(idon))./(parm.DONstd.^2*sum(dVt(iwet))));
Wo = Wo*scale;

px = zeros(nwet,11);
px(:,1:4) = Px(1:nwet,:);
pos = parm.pos;
pxx = zeros(nwet,pos(end,end));
posP = [1 2 3 4; 0 5 6 7; 0 0 8 9; 0 0 0 10];
for i1 = 1:4
    for i2 = i1:4
        pxx(:,pos(i1,i2)) = Pxx(1:nwet,posP(i1,i2));
    end
end
Ox = Nx(2*nwet+1:end,:); Oxx = Nxx(2*nwet+1:end,:);
P = DIP(iwet(ipo4)); px = px(ipo4,:); pxx = pxx(ipo4,:);
N = DIN(iwet(ino3)); nx = Nx(ino3,:); nxx = Nxx(ino3,:);
O = DON(iwet(idon)); ox = Ox(idon,:); oxx = Oxx(idon,:);

ep = P-PO4;   % DIP error beteen mod. and obs.;
en = N-NO3;   % DIN error beteen mod. and obs.;
eo = O-DON_o; % DON error beteen mod. and obs.;

f = 0.5*(ep.'*Wp*ep + en.'*Wn*en + eo.'*Wo*eo); % objective function;

fprintf('objective function value f = %3.3e \n',f);

% first derivative towards parameters;
if (nargout>1)
    fx = ep.'*Wp*px+en.'*Wn*nx + eo.'*Wo*ox;
end

% second derivatives towards parameters;
if (nargout>2)
    fxx = zeros(11,11);
    for i1 = 1:11
        for i2 = i1:11
            fxx(i1,i2) = px(:,i1).'*Wp*px(:,i2)+...
                ep.'*Wp*pxx(:,pos(i1,i2)) +...
                nx(:,i1).'*Wn*nx(:,i2)+...
                en.'*Wn*nxx(:,pos(i1,i2))+...
                ox(:,i1).'*Wo*ox(:,i2)+...
                eo.'*Wo*oxx(:,pos(i1,i2));
            fxx(i2,i1) = fxx(i1,i2); % symmetric
        end
    end    
end
toc
