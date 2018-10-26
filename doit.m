global GN
load po4obs_90x180x24.mat po4obs    % WOA2013 interpreted DIP
load no3obs_90x180x24.mat no3obs    % WOA2013 interpreted DIN 
po4int = po4obs; no3int = no3obs;   
clear po4obs no3obs
load raw_no3obs_90x180x24.mat no3raw % WOA2013 un-interpreted DIN
load raw_po4obs_90x180x24.mat po4raw % WOA2013 un-interpreted DIP
load o2obs_90x180x24 o2obs % WOA2013 un-interpreted O2
load npp_90x180.mat npp    % CbPM NPP
load transport_v4.mat      % OCIM transport operator
load DIN_dep_90x180.mat    % Atmospheric N deposition
load River_NP_Atm_P.mat output % Riverine N inputs
load DOMobs_wRefr_90x180x24.mat % DOM observation data.
%load N_OPT_N2P_slope.mat
load PN_DON1e-01_SDN100_preind.mat % a model N P field for initial
                                   % guess of Newton Method
DIN = real(DIN);
DON = real(DON);
PON = real(PON);
% get rid of negative values
ineg = find(no3raw(:)<0);
no3raw(ineg) = 0.05;
ineg = find(po4raw(:)<0);
po4raw(ineg) = 0.001;

% prescribe some constants.
dpa      = 365;       % days per year
spd      = 24*60^2;   % seconds per day
spa      = dpa*spd;   % seconds per year
n_iocn   = length(iocn); % number of wet gridboxes
I        = speye(n_iocn);
TRdiv    = -TR;     % [s^-1].
grd      = grid;

parm.TRdiv = TRdiv;
parm.M3d = M3d;
parm.grd = grd;

%%%%%%%%%%%% get dep data done %%%%%%%
SRF = M3d(:,:,1);
isrf = find(SRF(:));

%get river data and convert units from nmol/cm2/s to mmol/m3/s
DIN_RIV_FLUX = output.DIN_RIV_FLUX/grd.dzt(1)*1e-2;
DON_RIV_FLUX = output.DON_RIV_FLUX/grd.dzt(1)*1e-2;
DONr_RIV_FLUX = output.DONr_RIV_FLUX/grd.dzt(1)*1e-2;
% get rid of bad data
DIN_RIV_FLUX(isnan(DIN_RIV_FLUX))   = 0;
DON_RIV_FLUX(isnan(DON_RIV_FLUX))   = 0;
DONr_RIV_FLUX(isnan(DONr_RIV_FLUX)) = 0;
% check total amount
tDIN_RIV_tmp = DIN_RIV_FLUX.*dVt(:,:,1);
tDON_RIV_tmp = DON_RIV_FLUX.*dVt(:,:,1);
tDONr_RIV_tmp = DONr_RIV_FLUX.*dVt(:,:,1);
sum_tDIN_RIV = nansum(tDIN_RIV_tmp(:))*spa*14*1e-15;
sum_tDON_RIV = nansum(tDON_RIV_tmp(:))*spa*14*1e-15;
sum_tDONr_RIV = nansum(tDONr_RIV_tmp(:))*spa*14*1e-15;

[parm.DIN_RIV,parm.DON_RIV,parm.DONr_RIV] = deal(M3d*0);
parm.DIN_RIV(:,:,1) = 0.5*DIN_RIV_FLUX;
parm.DON_RIV(:,:,1) = 0.5*DON_RIV_FLUX;
parm.DIN_RIV = parm.DIN_RIV(iocn);
parm.DON_RIV = parm.DON_RIV(iocn);
parm.DONr_RIV(:,:,1) = DONr_RIV_FLUX;
fprintf(['riverine DIN, DON, DONr inputs are %3.3e, %3.3e, %3.3e ' ...
         'Tg/year,\n'], ...
        sum_tDIN_RIV,sum_tDON_RIV,sum_tDONr_RIV)

DIN_dep = NHy_pre_ind+NOx_pre_ind;   % nmol N/cm2/s
% DIN_dep = NHy_present+NOx_present;     % nmol N/cm2/s
DIN_dep = DIN_dep*1e-6*1e4/grd.dzt(1); % mmol N/m3/s
ibad = find(isnan(DIN_dep(isrf)) | isinf(DIN_dep(isrf)));
DIN_dep(isrf(ibad)) = 0;
% check total amount
tDN_tmp = DIN_dep.*dVt(:,:,1);
sum_tDN = nansum(tDN_tmp(:))*spa*14*1e-15;
parm.DIN_dep = M3d*0;
parm.DIN_dep(:,:,1) = DIN_dep;
parm.DIN_dep = parm.DIN_dep(iocn);
fprintf('DIN deposition is %3.3e Tg/year,\n',sum_tDN)

% DIP deposition is calculated based on dust.
%get river data and convert units from nmol/cm2/s to mmol/m3/s
DIP_RIV_FLUX = output.DIP_RIV_FLUX/grd.dzt(1)*1e-2;
DOP_RIV_FLUX = output.DOP_RIV_FLUX/grd.dzt(1)*1e-2;
DOPr_RIV_FLUX = output.DOPr_RIV_FLUX/grd.dzt(1)*1e-2;
DIP_RIV_FLUX(isnan(DIP_RIV_FLUX))   = 0;
DOP_RIV_FLUX(isnan(DOP_RIV_FLUX))   = 0;
DOPr_RIV_FLUX(isnan(DOPr_RIV_FLUX)) = 0;
% check total amount
tDIP_RIV_tmp = DIP_RIV_FLUX.*dVt(:,:,1);
tDOP_RIV_tmp = DOP_RIV_FLUX.*dVt(:,:,1);
tDOPr_RIV_tmp = DOPr_RIV_FLUX.*dVt(:,:,1);
sum_tDIP_RIV = nansum(tDIP_RIV_tmp(:))*spa*31*1e-15;
sum_tDOP_RIV = nansum(tDOP_RIV_tmp(:))*spa*31*1e-15;
sum_tDOPr_RIV = nansum(tDOPr_RIV_tmp(:))*spa*31*1e-15;

[parm.DIP_RIV,parm.DOP_RIV,parm.DOPr_RIV] = deal(M3d*0);
parm.DIP_RIV(:,:,1) = DIP_RIV_FLUX;
parm.DOP_RIV(:,:,1) = DOP_RIV_FLUX;
parm.DOPr_RIV(:,:,1) = DOPr_RIV_FLUX;
fprintf(['riverine DIP, DOP, DOPr inputs are %3.3e, %3.3e, %3.3e ' ...
         'Tg/year,\n'], sum_tDIP_RIV,sum_tDOP_RIV,sum_tDOPr_RIV)


Dust = output.dust_FLUX_IN; % g D/cm2/s
Dust = Dust*(1.0/0.98); % account for 2% dust lost
Dust = Dust*1e4; % g D/m2/s;
                 % from g D/m2/s to molP/m3/s
PO4_dep = Dust*0.00105*0.15/30.974/grd.dzt(1);
ibad = find(isnan(PO4_dep(isrf)) | isinf(PO4_dep(isrf)));
PO4_dep(isrf(ibad)) = 0;
% check total amount 96.5 Gg P yr^−1 (Mahowald et al 2008)
tDP_tmp = PO4_dep.*dVt(:,:,1);
sum_tDP = nansum(tDP_tmp(:))*spa*30.974*1e-9;
parm.DIP_dep = M3d*0;
parm.DIP_dep(:,:,1) = PO4_dep;
fprintf('DIP deposition is %3.3e Gg/year,\n',sum_tDP)
%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%

parm.c2n      = 106/16;  % C:N ratio
parm.p2c      = 0.006+0.0069*po4int; % P:C ratio from GB15 used to
                                     % scale NPP
parm.c2p      = 1./parm.p2c;

parm.Nc       = 1.8;          % critical DIN.
parm.sigma    = 0.30;
parm.kappa_p  = 1/(720*60^2); % DON remi.
parm.kappa_g  = 1/(1e6*spa);  % geological restoring.

ipo4 = find(~isnan(po4raw(iocn))); % index for valid DIP measurements
ino3 = find(~isnan(no3raw(iocn))); % index for valid DIN measurements

parm.o2obs    = o2obs;   
parm.POC      = PON(iocn)*parm.c2n;
parm.DIPbar   = sum(po4int(iocn).*dVt(iocn))./sum(dVt(iocn));
parm.DIPstd   = std(po4int(iocn));
parm.DINbar   = sum(no3int(iocn).*dVt(iocn))./sum(dVt(iocn));
parm.DINstd   = std(no3int(iocn));
parm.po4raw   = po4raw;
parm.no3raw   = no3raw;
parm.no3int   = no3int;
parm.po4obs   = po4int;
parm.DONobs   = DONobs;

parm.dVt = dVt;
% GN: global variable to reset the initial N-cycle model state for the
% Newton solver
GN  = 0.999*[DIN(iocn);PON(iocn);DON(iocn)];                                            
%%%%%%% prepare NPP for the model %%%%%%%%
inan = find(isnan(npp(:)));
npp(inan) = 0;
parm.npp = npp/(12*spd);
parm.Lambda = M3d*0;
parm.Lambda(:,:,1) = 0.5*(1/grd.dzt(1))*(1./parm.c2p(:,:,1))./(1e-9+po4int(:,:,1));
parm.Lambda(:,:,2) = 0.5*(1/grd.dzt(2))*(1./parm.c2p(:,:,2))./(1e-9+po4int(:,:,2));
parm.Lambda(:,:,3:end) = 0;
parm.Oc = 0;
parm.sl = 0.19;
parm.iN = 0.06;
parm.cc = 1; % parameter that scales water column oxygen.
parm.DIPc = 0;
%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load xhat_DON1e-01_SDN100_preind.mat % inital parameter guess
x0 = xhat;
nip = length(x0); % total number of adjustable parameters
% order of second derivatives
in = 1:11;
nin = length(in);
n = (nin+1)*nin/2;
pos = tril(ones(nin),0);
pos(pos==1) = 1:n;
pos = pos';
parm.in = in;
parm.nin = nin;
parm.pos = pos;

SDN_scale = 1.0; % scaling factor for sedimentary denitrification
                 % corresponding to "s" in the paper.

for jj = 1:length(SDN_scale)
    
    parm.SDN_scale = SDN_scale(jj);
    parm.DON_scale = 0.1;   
    [ib,o2] = buildD(parm,M3d);
    parm.ib = ib;
    parm.o2 = o2;
    myfun = @(x) neglogpost(x,parm);
    options = optimoptions(@fminunc,...
                           'Algorithm','trust-region',...
                           'GradObj','on',...
                           'Hessian','on',...
                           'Display','iter',...
                           'MaxFunEvals',2000, ...
                           'MaxIter',2000,...
                           'TolX',1e-9,...
                           'TolFun',1e-9,...
                           'DerivativeCheck','off',...
                           'FinDiffType','central',...
                           'PrecondBandWidth',Inf);
    
    [xhat,fval,exitflag] = fminunc(myfun,x0,options);
    
    [f,fx,fxx] = neglogpost(xhat,parm);
    fname = sprintf('xhat_DON%1.0e_SDN%d_preind', parm.DON_scale, parm.SDN_scale*100);
    save(fname,'xhat','f','fx','fxx')
    
    x0 = xhat;
    
end

%%%%%%%% test 1st and 2nd Derivative. %%%%%%%
test_the_code = 0;
if(test_the_code)
    dx = sqrt(-1)*eps.^3*eye(nip);
    for ii = 1:nip
        x  = real(x0)+dx(:,ii);
        parm.bp       = exp(x(1));
        parm.kappa_dp = exp(x(2));
        parm.alpha    = exp(x(3));
        parm.beta     = exp(x(4));
        parm.kw       = exp(x(5));
        parm.bn       = exp(x(6));
        parm.kappa_dn = exp(x(7));
        parm.cc       = exp(x(8));
        parm.n2p_l    = exp(x(9));
        parm.S        = exp(x(10));
        parm.Nc       = x(end);
        
        [f,fx,fxx] = neglogpost(x,parm);
        fprintf('%i %e %e %e %e %e %e %e %e %e %e %e %e \n',ii,real(fxx(ii,:))-imag(fx(:))'/eps^3);
    end
end