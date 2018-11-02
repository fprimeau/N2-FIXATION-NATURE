function [R,Rn,rhsRx,Rx,rhsRxx] = fixit(N,parm,Nx)
Nc = parm.Nc;
iwet = find(parm.M3d(:));
n_wet = length(iwet);
DIN = N(1:n_wet);

T = tanh(DIN-parm.Nc);
R = 0.5+0.5*T;         % corresponding to Eq. 6 

% calculate first and second derivatives;
if(nargout>1)
    dTdDIN = 1-T.^2;
    Rn = 0.5*dTdDIN;
end
if(nargout>2)
    dTdNc = -(1-T.^2);
    rhsRx = 0.5*dTdNc;
end
if(nargout>3)
    DINx = Nx(1:n_wet,:);
    Rx = d0(Rn)*DINx;
    Rx(:,11) = Rx(:,11)+rhsRx;

    % order of second derivatives
    nin = parm.nin;
    pos = parm.pos;

    % include all terms except those that involve Nxx
    dTdNc  = -(1-T.^2);
    dRdNc  = 0.5*dTdNc;

    % dNc dDIN
    d2TdNcdDIN = 2*T.*dTdDIN;
    d2RdNcdDIN = 0.5*d2TdNcdDIN;
    % dNc dNc
    d2TdNc2 = 2*T.*(T.^2 - 1);
    d2RdNc2 = 0.5*d2TdNc2;
    % dDIN dDIN
    d2TdDIN2 = 2*T.*(T.^2-1);
    d2RdDIN2 = 0.5*d2TdDIN2;
    %
    rhsRxx = zeros(n_wet,pos(end,end));
    rhsRxx(:,pos(11,11)) = d2RdNc2+d2RdNcdDIN.*DINx(:,11);
    % dNc dx
    for is = 1:11
        rhsRxx(:,pos(is,11)) = rhsRxx(:,pos(is,11))+d2RdNcdDIN.*DINx(:,is);
    end
    % dx dx
    for i1 = 1:nin
        for i2 = i1:nin
            rhsRxx(:,pos(i1,i2)) = rhsRxx(:,pos(i1,i2))+...
                DINx(:,i1).*d2RdDIN2.*DINx(:,i2);
        end
    end
end