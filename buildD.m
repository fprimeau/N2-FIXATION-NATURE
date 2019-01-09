function [ib,o2obs_c] = buildD(p,M3d)
    
% the following function is used to build denitrification operator
% by assuming global denitrification amount is fixed at 150 TgN/year
% and water colume denitrification is one third of total denitrification
% we also assume denitrification is propotional to oxygen deficit
% and POC concentration.
    o2obs = p.o2obs;
    iocn = find(M3d(:));
    n_iocn = length(iocn);
    I = speye(n_iocn);
    [nx,ny,nz] = size(M3d);
    msk = M3d(:,:,[2:end,end]);
    B = M3d-msk;
    B(:,:,end) = M3d(:,:,24);       % bottom layer.

    %ib = B(iocn);                   % bottom layer for the whole ocean.
    B(:,:,1:6) = 0;                 % get rid of surface bottom grid;
    iB = B(iocn);                   % bottom layer for deep ocean (>300m)
    
    S = M3d*0;
    junk1 = M3d;
    junk1(:,:,8:end) = 0;             % get upper ocean (300 m) grid mask
    junk2 = junk1(:,:,[2:end,end]);   % shift grid box up.
    S = junk1-junk2;                 % get bottom layer
    S(:,:,7) = 0;                    
    iS = S(iocn);                    % bottom layer for surface ocean.
    
    ib = (iB + iS) * p.SDN_scale;
    % ib =  iB + iS  * p.SDN_scale;
    o2obs_c = M3d*0;
    o2obs(iocn) = o2obs(iocn).*44.661;        % convert unit form [ml/l] to [umol/l].
    o2obs_c(iocn) = o2obs(iocn).*1.009-2.523; % o2 correction based on Bianchi et al.(2012) [umol/l] .
    ineg = find(o2obs_c(:)<0);                % find out negative values and set them to zero.
    o2obs_c(ineg) = 0;
    

