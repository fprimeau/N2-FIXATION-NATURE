function [PFdiv,dPFDdb,d2PFDdb2,dPFDdr,d2PFDdr2,d2PFDdbr] = buildPFD(M3d,b,r,grd);
% this code is used to build Particle Flux Diverngence (PFD)
%
%                      ________________________                                        
%                     /         A             /|  POP sinking          
%                 top/_______________________/ |       |       
%                    |  |                    | |       | -w   
%                    | /        V            | /       |       
%                 bot|/______________________|/        V
%                                                  
% PFD = (A*w(top)*POP(top)-A*w(bot)*POP(bot))/V;
% add an exra layer of zeros at the bottom to ensure there is no
% flux in or out of the domain when using periodic shift operators
% for finite differences and averaging
    [ny,nx,nz] = size(M3d);
    M3D = zeros(ny,nx,nz+1);
    M3D(:,:,1:end-1) = M3d;
    % add the zw coordinate at the top of the extra layer
    ZW3d = grd.ZW3d;
    ZW3d = ZW3d(:,:,[1:end,end]);
    ZW3d(:,:,end) = grd.ZW3d(:,:,end)+grd.dzt(end);
    % areas of the top of the grid box
    dAt = (grd.DXT3d.*grd.DYT3d).*M3d;
    % volume of the grid boxes
    dVt = (grd.DXT3d.*grd.DYT3d.*grd.DZT3d).*M3d;
    %
    n = nx*ny*(nz+1);
    I0 = speye(n);
    i0 = zeros(ny,nx,nz+1);
    i0(:) = 1:n;
    % periodic shifts OK because M3D has a layer of zeros on the bottom
    iu = i0(:,:,[nz+1,1:nz]); %use a periodic upward shift
    ib = i0(:,:,[2:nz+1,1]); % use a periodic downward shift
    IU = I0(iu,:);
    IB = I0(ib,:);
    % keep only wet boxes
    iocn = find(M3D(:));
    I0 = I0(iocn,:); I0 = I0(:,iocn);
    IU = IU(:,iocn); IU = IU(iocn,:);
    IB = IB(:,iocn); IB = IB(iocn,:);
    % (averages POP onto the top of the grid boxes)
    AVG = d0((I0+IU)*M3D(iocn))\(I0+IU);
    % (compute the divergence in the center of the grid boxes)
    DIV = d0(dVt(iocn))\(I0-IB)*d0(dAt(iocn));
    % (compute the flux at the top of the grid cells)
    % mimics a Martin curve flux attenuation profile
    %(see Kriest and Oschelies 2008 in Biogeosciences)
    a = r/b;
    % particle sinking velocity at the top of the grid cells.
    MSK = M3D.*M3D(:,:,[nz+1,1:nz]);
    M = MSK.*ZW3d;
    w = -a*M;
    dadb = -r/(b^2);
    dadr = 1/b;
    dwdb = -dadb*M;
    dwdr = -dadr*M;
    d2adb2 = 2*r/(b^3);
    d2adr2 = 0;
    d2adbr = -1/b^2;
    d2wdb2 = -d2adb2*M;
    d2wdr2 = -d2adr2*M;
    d2wdbr = -d2adbr*M;
    %FLUX = d0(w(iocn))*AVG;
    FLUX = d0(w(iocn))*IU;
    dFLUXdb = d0(dwdb(iocn))*IU;
    dFLUXdr = d0(dwdr(iocn))*IU;
    d2FLUXdb2 = d0(d2wdb2(iocn))*IU;
    d2FLUXdbr = d0(d2wdbr(iocn))*IU;
    d2FLUXdr2 = 0;
    % particle flux divergence operator
    PFdiv = DIV*FLUX;
    dPFDdb = DIV*dFLUXdb;
    dPFDdr = DIV*dFLUXdr;
    d2PFDdb2 = DIV*d2FLUXdb2;
    d2PFDdbr = DIV*d2FLUXdbr;
    d2PFDdr2 = DIV*d2FLUXdr2;
end 
