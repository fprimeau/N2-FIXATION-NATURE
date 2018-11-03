global GN

% OCIM transport operator including grid (grid) and mask (M3d) and
% advection and diffusion transport operator (TR);
load transport_v4.mat   
load envi_data.mat no3obs po4obs o2obs npp % environmental data

% raw uninterpreted data used to constrain the model
load constraint_data po4raw no3raw DONobs   

% 2D riverine and atmospheric sources based on Lamarque et al.,
% 2010 and Seitzinger et al 2010.
load external_sources.mat  

% load N P fields for initial iterate of Newton Method
load initial_iterate.mat DIN DON PON
grd   = grid;         % grid info. for OCIM model;
iwet  = find(M3d(:)); % find index of wet grid box;
nwet  = length(iwet); % number of wet grid boxes;
I     = speye(nwet);  % build an identity matrix;

% global vector to pass initial N guess to N function.
GN  = 0.999*[DIN(iwet);PON(iwet);DON(iwet)];

% get rid of negative values
ineg = find(no3raw(:)<0);
no3raw(ineg) = 0.05;
ineg = find(po4raw(:)<0);
po4raw(ineg) = 0.01;

% prescribe some constants.
dpa = 365;      % days per year;
spd = 24*60^2;  % seconds per day;
spa = dpa*spd;  % seconds per year;

% OCIM model info. loaded from transport_v4;
parm.M3d = M3d;   % land ocean mask;
parm.TRdiv = -TR; % advection/diffusion transport operator [s^-1];
parm.grd = grd;   

%%%%%%%%%%%% get external source data done %%%%%%%
% find surface wet box index;
SRF = M3d(:,:,1);     % get surface layer
isrf = find(SRF(:));  % find indices for surface wet grid boxes;

%get river data and convert units from nmol/cm2/s to mmol/m3/s
DIN_RIV_FLUX = DIN_RIV_FLUX/grd.dzt(1)*1e-2;
DON_RIV_FLUX = DON_RIV_FLUX/grd.dzt(1)*1e-2;
DONr_RIV_FLUX= DONr_RIV_FLUX/grd.dzt(1)*1e-2;

% get rid of bad data
DIN_RIV_FLUX(isnan(DIN_RIV_FLUX))   = 0;
DON_RIV_FLUX(isnan(DON_RIV_FLUX))   = 0;
DONr_RIV_FLUX(isnan(DONr_RIV_FLUX)) = 0;

% check and print total amount, where dVt is a 3D gridbox volume matrix;
tDIN_RIV_tmp = DIN_RIV_FLUX.*dVt(:,:,1);
tDON_RIV_tmp = DON_RIV_FLUX.*dVt(:,:,1);
tDONr_RIV_tmp = DONr_RIV_FLUX.*dVt(:,:,1);
sum_tDIN_RIV = nansum(tDIN_RIV_tmp(:))*spa*14*1e-15;
sum_tDON_RIV = nansum(tDON_RIV_tmp(:))*spa*14*1e-15;
sum_tDONr_RIV = nansum(tDONr_RIV_tmp(:))*spa*14*1e-15;
fprintf(['riverine DIN, DON, DONr inputs are %3.3e, %3.3e, %3.3e ' ...
         'Tg/year,\n'], ...
        sum_tDIN_RIV,sum_tDON_RIV,sum_tDONr_RIV)

% assign external data to 3D grid;
[parm.DIN_RIV,parm.DON_RIV,parm.DONr_RIV] = deal(M3d*0);
parm.DIN_RIV(:,:,1) = 0.5*DIN_RIV_FLUX; % assume 50% riverine inputs
parm.DON_RIV(:,:,1) = 0.5*DON_RIV_FLUX; % reach open ocean;
parm.DIN_RIV = parm.DIN_RIV(iwet);
parm.DON_RIV = parm.DON_RIV(iwet);
parm.DONr_RIV(:,:,1) = DONr_RIV_FLUX;

% get atmospheric N deposition
DIN_dep = NHy_pre_ind+NOx_pre_ind;   % nmol N/cm2/s, preindustrial;
% DIN_dep = NHy_present+NOx_present;    % nmol N/cm2/s, modern day;
DIN_dep = DIN_dep*1e-6*1e4/grd.dzt(1); % mmol N/m3/s, unit convertion;
ibad = find(isnan(DIN_dep(isrf)) | isinf(DIN_dep(isrf))); % remove
                                                          % bad date
DIN_dep(isrf(ibad)) = 0;

% check total amount
tDN_tmp = DIN_dep.*dVt(:,:,1);
sum_tDN = nansum(tDN_tmp(:))*spa*14*1e-15;
parm.DIN_dep = M3d*0;
parm.DIN_dep(:,:,1) = DIN_dep;
parm.DIN_dep = parm.DIN_dep(iwet);
fprintf('DIN deposition is %3.3e Tg/year,\n',sum_tDN)

% DIP deposition is calculated based on dust.
%get river data and convert units from nmol/cm2/s to mmol/m3/s
DIP_RIV_FLUX  = DIP_RIV_FLUX/grd.dzt(1)*1e-2;
DOP_RIV_FLUX  = DOP_RIV_FLUX/grd.dzt(1)*1e-2;
DOPr_RIV_FLUX = DOPr_RIV_FLUX/grd.dzt(1)*1e-2;
DIP_RIV_FLUX(isnan(DIP_RIV_FLUX))   = 0;
DOP_RIV_FLUX(isnan(DOP_RIV_FLUX))   = 0;
DOPr_RIV_FLUX(isnan(DOPr_RIV_FLUX)) = 0;

% check and print total amount
tDIP_RIV_tmp = DIP_RIV_FLUX.*dVt(:,:,1);
tDOP_RIV_tmp = DOP_RIV_FLUX.*dVt(:,:,1);
tDOPr_RIV_tmp = DOPr_RIV_FLUX.*dVt(:,:,1);
sum_tDIP_RIV = nansum(tDIP_RIV_tmp(:))*spa*31*1e-15;
sum_tDOP_RIV = nansum(tDOP_RIV_tmp(:))*spa*31*1e-15;
sum_tDOPr_RIV = nansum(tDOPr_RIV_tmp(:))*spa*31*1e-15;
fprintf(['riverine DIP, DOP, DOPr inputs are %3.3e, %3.3e, %3.3e ' ...
         'Tg/year,\n'], sum_tDIP_RIV,sum_tDOP_RIV,sum_tDOPr_RIV)

% assign external data to 3D grid;
[parm.DIP_RIV,parm.DOP_RIV,parm.DOPr_RIV] = deal(M3d*0);
parm.DIP_RIV(:,:,1) = DIP_RIV_FLUX;
parm.DOP_RIV(:,:,1) = DOP_RIV_FLUX;
parm.DOPr_RIV(:,:,1) = DOPr_RIV_FLUX;

Dust = dust_FLUX_IN; % g D/cm2/s
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

parm.c2n = 106/16;              % C:N ratio
parm.p2c = 0.006+0.0069*po4obs; % P:C ratio (GB15) used to scale NPP

% wrap data into structure for easy passing to other functions;
parm.o2obs    = o2obs;   % oxygen observation used to scale WC denitrification;
parm.POC      = PON(iwet)*parm.c2n;
parm.DIPbar   = sum(po4obs(iwet).*dVt(iwet))./sum(dVt(iwet));
parm.DIPstd   = std(po4obs(iwet));
parm.DINbar   = sum(no3obs(iwet).*dVt(iwet))./sum(dVt(iwet));
parm.DINstd   = std(no3obs(iwet));
parm.po4raw   = po4raw;
parm.no3raw   = no3raw;
parm.no3obs   = no3obs;
parm.po4obs   = po4obs;
parm.DONobs   = DONobs;

parm.dVt = dVt;

%%%%%%% prepare NPP for the model %%%%%%%%
inan = find(isnan(npp(:)));
npp(inan) = 0;
parm.npp = npp/(12*spd);
parm.Lambda = M3d*0;
parm.Lambda(:,:,1) = 0.5*(1/grd.dzt(1))*parm.p2c(:,:,1)./(1e-9+po4obs(:,:,1));
parm.Lambda(:,:,2) = 0.5*(1/grd.dzt(2))*parm.p2c(:,:,2)./(1e-9+po4obs(:,:,2));
parm.Lambda(:,:,3:end) = 0;
%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% prescribed parameters %%%%%%%%%%%%%%
parm.sigma    = 0.30;         % fraction of NPP directed to DOP;
parm.kappa_p  = 1/(720*60^2); % DON remi;
parm.kappa_g  = 1/(1e6*spa);  % geological restoring;
parm.DIPc = 0;  % corresponding to DIPc in Eq. S7;
parm.Oc = 0;    % corresponding to x_c in Eq. S4;
parm.sl = 0.19; % parameter for Benthic denitrification (Bohlen et al.,2012);
parm.iN = 0.06; % parameter for Benthic denitrification (Bohlen et al.,2012);
%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% parameters need to be optimized %%%%%%%
bp       = 0.96;     % corresponds to b_P in table S1;
kappa_dp = 0.94e-7;  % corresponds to kappa_dP in table S1;
alpha    = 1.71e-2;  % corresponds to alpha in table S1;
beta     = 0.57;     % corresponds to beta in table S1;
kw       = 0.91e-9;  % corresponds to k_w in table S1;
bn       = 0.80;     % corresponds to b_N in table S1;
kappa_dn = 0.57e-7;  % corresponds to kappa_dN in table S1;
cc       = 20.64;    % corresponds to Delta in table S1;
n2p_l    = 12.41;    % corresponds to A in table S1;
S        = 6.58;     % corresponds to B in table S1;
Nc       = -1.23;    % corresponds to [DIN]_c in table S1;
%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%

% log transform parameters that have to be positive;
x0 = [log([bp; kappa_dp; alpha; beta; kw; bn; ...
           kappa_dn; cc; n2p_l; S]); Nc];

nip = length(x0); % number of optimized parameters;

% order of second derivatives
in = 1:nip; % parameter indices
nin = length(in);
n = (nin+1)*nin/2;       % number of unique 2nd derivatives
                         
% the pos matrix is an nip x nip array with numbers from
% 1 to n in the upper triangular part. These numbers are used as
% the column index in a big array whos columns store the 2nd
% derivative of the model state with respect to the parameters
pos = tril(ones(nin),0); 
pos(pos==1) = 1:n;
pos = pos'; 
parm.in = in;
parm.nin = nin;
parm.pos = pos;

SDN_scale = 1.0; % scaling factor for sedimentary denitrification
                 % corresponds to "s" in the paper.

for jj = 1:length(SDN_scale)
    
    parm.SDN_scale = SDN_scale(jj); 
    parm.DON_scale = 0.1;           % parameter that scales DON weight;
    [ib,o2] = buildD(parm,M3d);     % find bottom grid box and
                                    % correct O2 concentration
                                    % based on Bianchi et al., (2012);
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
%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%