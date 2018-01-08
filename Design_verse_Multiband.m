clear all;
close all;

AM_only = 0;
Nt = 2048;
gamma_mT = 2*pi*4.257*1e4; %<--- same as in minTime gradient function
slthick = 2*1e-3;
b1max = 13*1e-3; % <--- peak B1 [mT]
maxg = 40; %<-- Maximum gradient amplitude [mT/m]
maxgslew = 200*1e3; % <--- Maximum gradient amplitude [mT/m/s]

% --- Pulse characteristics --- %
mb = 4;
tb = 4;
bs = 14;

% --- Select GIRF --- %
girf_idc = 1;

switch girf_idc
    case 1
    % h1_GIRF_20140729 is a measured GIRF used for experimental results.
    girf = load('bin/h1_GIRF_20140729');

    case 2
    % h2_GIRF_20170901 is a reconstructed GIRF from Testud, using an
    % alternative vendor gradient system with a higher temporal BW.
    girf = load('bin/h2_GIRF_20170901');
    
    case 3

    % Thirdly, it's possible to specify an analytical GIRF model. Here, a
    % mono-exponential filter will roughly approximate h1 to have the frequency
    % response of as a Lorentzian, with a fixed delay-term.

    tau = 41.717*1e-6; %<-- Time-constant for mono-exponential filter
    delay =  45.442*1e-6; %<-- constant delay.
    girf = [tau delay];
end

%% Design Singleband pulse
load('SB_SLR_cvxdesign_flip180_quad_Mar27.mat')
rfsb = pulse(tb-1).rf;                    

rfsb = length(rfsb)/Nt*interp1(linspace(0,1,length(rfsb)),rfsb,linspace(0,1,Nt))';                    

dt_sb = max(abs(rfsb))/(gamma_mT*b1max);
rfsb_mT = rfsb ./(gamma_mT*dt_sb);

BW_sb = tb/(length(rfsb_mT)*dt_sb);
Gsel_sb = 2*pi*BW_sb/(gamma_mT*slthick);                                      
Gz_sb = Gsel_sb*ones(length(rfsb_mT),1);

%% Design Phase-optimized MB pulse
rfmb = Phaseopt_fn(rfsb,mb,tb,bs,AM_only);

dt=max(abs(rfmb))./(gamma_mT*b1max);
% Convert MB pulse from rad to mT:
rfmb = rfmb./(gamma_mT*dt);
BW = tb/(Nt*dt);
Gsel = 2*pi*BW/(gamma_mT*slthick);
Gz = Gsel*ones(Nt,1);
gmb = [0*Gz 0*Gz Gz];    
t = (0:length(rfmb)-1)*dt;


%% Design linear-phase MBv pulse

fprintf('Designing MBv pulse...\n')
verse_singleband = 0;
dt_os = 2;

[rfMBv,gMBv,gMBv_actual]= dz_MBverseGirf(rfmb,Gz,dt,maxg,...
    maxgslew,b1max,verse_singleband,mb,bs*slthick,dt_os,AM_only,girf);
Nv = length(rfMBv);
dtv =dt/dt_os;
tv = 0:dtv:(Nv-1)*dtv;

%% Design linear-phase vMB pulse
fprintf('Designing vMB pulse...\n')

rf_init = rfsb_mT;

% Scale input single-band RF to mT:
verse_singleband = 1;
dt_os = 2;
epsilon = 1e-4;
max_iterations = 1;

[rfvMB,gvMB,gvMB_actual]= dz_MBverseGirf(rf_init,Gz_sb,dt_sb,maxg,...
    maxgslew,b1max,verse_singleband,mb,bs*slthick,dt_os,AM_only,girf);
Nvmb = length(rfvMB);
dtvmb =dt_sb/dt_os;
tvmb = 0:dtvmb:(Nvmb-1)*dtvmb;

%% Run Bloch simulations.
fprintf('\nRunning Bloch simulations...\n')
spos = (1:mb)-(mb+1)/2; 

Nz = 4096;
FOV = 3*(floor(mb/2)+1)*bs*slthick;
z = linspace(-FOV/2,FOV/2,Nz)';

pos = [z(:)*0 z(:)*0 z(:)];

mxy_fa =@(a,b)180/pi*acos(a.*conj(a)-b.*conj(b));
mxy_ref = @(a,b)b.^2;

% --- Select whether to plot refocusing profiles or in pure Flip-angle ---
% %
plot_refocusing = 0;
if plot_refocusing
    mxy_display = mxy_ref;
else
    mxy_display = mxy_fa;
end

% Simulate Constant gradient MB pulse
fprintf('Running for constant gradient MB pulses...\n')
[~,~,~,~,a,b] = blochsim_CK(rfmb,gmb,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dt);
mxy0 = mxy_display(a(:,end),b(:,end));

fprintf('Running for MBv pulses ...\n')
[~,~,~,~,a,b] = blochsim_CK(rfMBv,gMBv,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dtv);
mxy_mbv = mxy_display(a(:,end),b(:,end));
[~,~,~,~,a,b] = blochsim_CK(rfMBv,gMBv_actual,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dtv);
mxy_mbv_act = mxy_display(a(:,end),b(:,end));

fprintf('Running for vMB pulses...\n')
[~,~,~,~,a,b] = blochsim_CK(rfvMB,gvMB,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dtvmb);
mxy_vmb = mxy_display(a(:,end),b(:,end));
[~,~,~,~,a,b] = blochsim_CK(rfvMB,gvMB_actual,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dtvmb);
mxy_vmb_act = mxy_display(a(:,end),b(:,end));

%% Plot RF, gradient waveforms and slice-profiles
fh = figure;
nr = 3;
nc = 5;

% Plot constant-gradient MB
subplot(nr,nc,0*nc+1);plot(t*1e3,abs(rfmb)*1e3);
ylabel('MB (const)');

subplot(nr,nc,0*nc+2);plot(t*1e3,abs(gmb(:,3)));

subplot(nr,nc,0*nc+3:nc);plot(z*1e2,abs(mxy0));
legend('Target');

% Plot MB verse
subplot(nr,nc,1*nc+1);plot(tv*1e3,abs(rfMBv)*1e3);
ylabel('MBv');

subplot(nr,nc,1*nc+2);plot(tv*1e3,abs(gMBv(:,3)));
hold on;plot(tv*1e3,abs(gMBv_actual(:,3)));

subplot(nr,nc,1*nc+3:2*nc);plot(z*1e2,abs(mxy_mbv));
hold on;plot(z*1e2,abs(mxy_mbv_act));
legend('Target','Predicted gradient distortion');

% Plot verse MB
subplot(nr,nc,2*nc+1);plot(tvmb*1e3,abs(rfvMB)*1e3);
ylabel('vMB');

subplot(nr,nc,2*nc+2);plot(tvmb*1e3,abs(gvMB(:,3)));
hold on;plot(tvmb*1e3,abs(gvMB_actual(:,3)));

subplot(nr,nc,2*nc+3:3*nc);plot(z*1e2,abs(mxy_vmb));
hold on;plot(z*1e2,abs(mxy_vmb_act));
legend('Target','Predicted gradient distortion');