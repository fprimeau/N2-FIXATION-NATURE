function [wcD,sdD,Dn,rhsDx,Dx,rhsDxx] = denit_don(N,PFdiv_n,dPFDdb,d2PFDdb2,parm,Nx)
% compute the rhs of the water column and sedimentary
% denitrification terms for the computation of Nxx

% unpack the parameters to be optimized
Oc = parm.Oc;
cc = parm.cc;
kw = parm.kw;
%
c2n = parm.c2n;
iwet = find(parm.M3d(:));
n_wet = length(iwet);
%
DIN = N(1:n_wet);
PON = N(n_wet+1:2*n_wet);
DON = N(2*n_wet+1:end);


% mask for where WC denitrification occurs
o2 = parm.o2(iwet);
[io,g1,g2,g11,g12,g22] = getio(o2,Oc,cc);

iB = parm.ib;
sl = parm.sl;
iN = parm.iN;
iX = sl*exp(log(0.99)*(o2-DIN))+iN;
diXdDIN = -sl*exp(log(0.99)*(o2-DIN))*log(0.99);
d2iXdDIN2 = sl*exp(log(0.99)*(o2-DIN))*log(0.99)^2;

wcD = c2n*kw*io.*DON.*DIN;
sdD = -c2n*(PFdiv_n*PON).*iB.*iX;

if(nargout>2)
    nin = parm.nin;
    Z = sparse(n_wet,n_wet);
    Dn = [d0(c2n*kw*io.*DON)-d0(c2n*(PFdiv_n*PON).*iB.*diXdDIN),...
          -d0(c2n*iB.*iX)*PFdiv_n,d0(c2n*kw*io.*DIN)];
          
    rhsDx = zeros(n_wet,nin);
    rhsDx(:,5) = c2n*kw*io.*DON.*DIN;
    rhsDx(:,6) = -c2n*(parm.bn*dPFDdb*PON).*iB.*iX;
    %rhsDx(:,8) = c2n*kw*g1.*DON.*DIN;
    rhsDx(:,8) = c2n*kw*g2.*DON.*DIN;

end    
if(nargout>4)
    DINx = Nx(1:n_wet,:);
    PONx = Nx(n_wet+1:2*n_wet,:);
    DONx = Nx(2*n_wet+1:end,:);

    % order of second derivatives
    pos = parm.pos;
    kwio = kw.*io;
    kwg = zeros(n_wet,pos(1,end));
    kwg(:,5) = kw*io;
    %kwg(:,8) = kw*g1;
    kwg(:,8) = kw*g2;

    kwgg = zeros(n_wet,pos(end,end));
    kwgg(:,pos(5,5)) = kw*io;
    %kwgg(:,pos(5,8)) = kw*g1;
    kwgg(:,pos(5,8)) = kw*g2;
    %kwgg(:,pos(8,8)) = kw*g11;
    %kwgg(:,pos(8,9)) = kw*g12;
    kwgg(:,pos(8,8)) = kw*g22;
    Wxx = zeros(n_wet,pos(end,end));
    for i1=1:11
        for i2 = i1:11
            Wxx(:,pos(i1,i2)) = Wxx(:,pos(i1,i2)) +...
                c2n*kwgg(:,pos(i1,i2)).*DON.*DIN +...
                c2n*kwg(:,i1).*DON.*DINx(:,i2)+...
                c2n*kwg(:,i2).*DON.*DINx(:,i1)+...
                c2n*kwg(:,i1).*DONx(:,i2).*DIN +...
                c2n*kwg(:,i2).*DONx(:,i1).*DIN+...
                c2n*kwio.*DONx(:,i1).*DINx(:,i2)+...
                c2n*kwio.*DONx(:,i2).*DINx(:,i1);
        end
    end
                
    Sxx = zeros(n_wet,pos(end,end));
    for i1=1:5
        Sxx(:,pos(i1,6)) = Sxx(:,pos(i1,6))+...
            d0(-c2n*iB.*iX)*(parm.bn*dPFDdb*PONx(:,i1))+...
            d0(-c2n*iB.*(parm.bn*dPFDdb*PON))*(d0(diXdDIN)*DINx(:,i1));
    end
    Sxx(:,pos(6,6)) = Sxx(:,pos(6,6))-...
          c2n*iB.*iX.*((parm.bn*dPFDdb+parm.bn^2*d2PFDdb2)*PON)-...
        2*c2n*iB.*iX.*(parm.bn*dPFDdb*PONx(:,6))-...
        2*c2n*iB.*(parm.bn*dPFDdb*PON).*diXdDIN.*DINx(:,6);
    for i1 = 7:nin
        Sxx(:,pos(6,i1)) = Sxx(:,pos(6,i1))+...
            d0(-c2n*iB.*iX)*(parm.bn*dPFDdb*PONx(:,i1))+...
            d0(-c2n*iB.*(parm.bn*dPFDdb*PON))*(d0(diXdDIN)*DINx(:,i1));
    end

    for i1=1:nin
        for i2 = i1:nin
            Sxx(:,pos(i1,i2)) = Sxx(:,pos(i1,i2))-...
                c2n*iB.*(PFdiv_n*PONx(:,i1)).*diXdDIN.*DINx(:,i2)-...
                c2n*iB.*(PFdiv_n*PONx(:,i2)).*diXdDIN.*DINx(:,i1);
        end
    end
    for is = 1:nin
        r = [pos(is,is):pos(is,end)];
            Sxx(:,r) = Sxx(:,r)-c2n*...
                d0(iB.*(PFdiv_n*PON).*d2iXdDIN2.*DINx(:,is))*DINx(:,is:end);
    end

    rhsDxx = Wxx+Sxx;

    wcDx = d0(c2n*kw*io.*DON)*DINx+...
           d0(c2n*kw*io.*DIN)*DONx;
    wcDx(:,5) = wcDx(:,5)+c2n*kw*io.*DON.*DIN;
    %wcDx(:,8) = wcDx(:,8)+c2n*kw*g1.*DON.*DIN;
    wcDx(:,8) = wcDx(:,8)+c2n*kw*g2.*DON.*DIN;
    
    sdDx = -d0(c2n*(PFdiv_n*PON).*iB.*diXdDIN)*DINx-...
           d0(c2n*iB.*iX)*(PFdiv_n*PONx);
    sdDx(:,6) = sdDx(:,6)-c2n*(dPFDdb*PON).*iB.*iX;

    Dx = wcDx+sdDx;

    
end
    
function [io,g1,g2,g11,g12,g22] = getio(o2,Oc,cc)
io = 0.5*(1-tanh((o2-Oc)/cc));
% syms o2 Oc cc
% g = 0.5*(1-tanh((o2-Oc)/cc))
% g1 = simplify(diff(g,Oc)*Oc);
% g2 = simplify(diff(g,cc)*cc);
% g11 = simplify(diff(g1,Oc)*Oc);
% g12 = simplify(diff(g1,cc)*cc);
% g21 = simplify(diff(g2,Oc)*Oc);
% g22 = simplify(diff(g2,cc)*cc);
g2 = 0.5*(sech((o2-Oc)/cc).^2).*(o2-Oc)/cc;
%
g1 = 0.5*sech((o2-Oc)/cc).^2*Oc/cc;
%
g12 = (Oc*(tanh((Oc - o2)/cc).^2 - 1).*(cc - 2*Oc*tanh((Oc - o2)/cc) + 2*o2.*tanh((Oc - o2)/cc)))/(2*cc^2);
%
g11 = (Oc*(cc - 2*Oc*tanh((Oc - o2)/cc)))./(2*cc^2*cosh((Oc - o2)/cc).^2);
%
g22 = -((Oc - o2).*(tanh((Oc - o2)/cc).^2 - 1).*(cc - 2*Oc*tanh((Oc - o2)/cc) + 2*o2.*tanh((Oc - o2)/cc)))/(2*cc^2);
