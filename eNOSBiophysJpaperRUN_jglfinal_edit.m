%% clear everything 
clc; close all; clear all;

%% global variables
global tau mu1 rr aktt pkct kcatp camtot bufftot pip2t phiinf kleak delta Kcoup SGC0 kpip2 kmpip2 kmakt kcam0 kmcam keakt kmeakt kpt5cam kon koff qmax
global kcce cs0 cex kres krel kout vr VpPKC kcavdec Gt ka kd alp api3k rstar W rcgmp Vdg kmdg NOdec RNO eta k2 k3 k5 k2p theta N xi kmpkc kdcam
global kca4cam a0 a1 g0 g1 b1 


%% Initialize parameters
tmax = 100000; % maximum time
tau = 0.0; % [dynes/cm^2] initialize 'basal' state of endothelial cells at 'zero' shear stress
% (putting a v small value like tau=10^-8 does not change initialization) 

mu1=0.2; rr=10.00;                             % degradation parameter and recycling parameter in IP3-PIP2 reactions
aktt=0.1; enost=0.04; pkct=aktt;               % total amounts in uM of enos, akt and pkc, assumed to be at ~ 10 nM concentrations, consistent with various expts cited
kcatp=0.026;                                   % in IP3 eqn from plank 2005 paper, atp term

camtot=30; bufftot=120; pip2t=10.0;    % Total amounts of Calmodulin, Caveolin, PIP2 and Buffer protein, in uM

phiinf=0.1;                                    % in IP3 eqn, the bulk ATP term in Plank's 2005 paper
kleak = 10^(-7);                               % Calcium leak term, discussed in paraemter initialization etc. 
delta=24;                                      % discussed in results section
lambdacap=15;                                  % parameter for shear-stress driven constitutive activation of GPCRs, discussed in appendix
Kcoup=0.002;                                   % Parameter to relate Lemon's model for GPCR activation, PIP2 etc and Plank's model. See appendix
SGC0=0.1;                                      % Soluble guanylate cyclase. See supplementary material

kpip2=0.021; kmpip2=0.7024; kmakt=0.1155;                                       % parameters for PIP eqns
kcam0=7.5; kmcam=0.01; keakt=0.004; kmeakt=keakt/18; kpt5cam=3;                 % parameters for akt and cam eqns, including for eNOS activation by both
kon=100; koff=300;                                                              % parameters for Ca buffering
qmax=17.6; kcce=8*10^-7; cs0=2828; cex=1500;                               % these are mostly parameters for the calcium accumulation eqns
kres=5; krel=6.64; kout=24.7; vr=3.5;                                           % parameters for calcium fluxes and resequestration
VpPKC=0.12/60; kcavdec=VpPKC*1/9;                      		% parameters for enos-pkc and cav eqns
Gt=0.332; %10^5;                                                                % Concentration of G-proteins in uM based on 10^5 molecules in a cell volume = 5*10^(-16) as per Lemon et al
ka = 0.017; kd=0.15;                                                            % kinetic activation parameters for g-proteins
alp=2.781/0.332;%*0.332/10^5;                                                   % alpha for IP3 activation; if using Gt = 10^5 molecules uncomment conversion factor

api3k=2.5;

rstar=tanh(3.14*tau/lambdacap);                                                 % shear stress activation of GPCRs
W=1.4*(tau/10+sqrt(16*2.86^2+tau^2/100)-4*2.86)^2/(tau/10+sqrt(16*2.86^2+tau^2/100));  % for calculation of plank et al's Ca2+ equations 

rcgmp=1.26; Vdg=0.0695; kmdg=2; NOdec=382; RNO=300.0;                          % parameters for NO production or consumption and cGMP production; refer appendix for RNO values

eta=0.0030000;                                                                 % 0.003 is the eta value used in the study as discussed in results section
k2=0.2; k3=0.15; k5=0.32; k2p=0.0022;
theta = .0045; %kca4cam equation
N = 2;
xi=.15/20;
kmpkc = kmakt;
kdcam = 1; %kca4cam equation
kca4cam = kon; %kca4cam equation (equal to kon)
a0 = 1200.1667*10^-6;
a1 = 37.34*10^-3;
g0 = 4.8006*10^-6;
g1 = 35.3*10^-3;
b1 = 15.15*10^-3;

prec=1e-10; % tolerances
options = odeset('RelTol',prec,'AbsTol',prec.*ones(16,1));

% calculates the 'basal' state of the EC at zero shear using matlab ode 
% equations are fairly stiff so we use ode15s()
[T,X0] = ode15s(@eNOSBiophysJpaperequations_jglfinal_edit,[0 tmax],[0, 0, pip2t, .0042, .1, cs0, 0, .0031, .0031, 0,0,0,0,0,0,enost],options);   
initValues = X0(end,:);

%% tau == 8 dyn/cm^2
tau=8;   % non-zero value of tau, we now check how much a non-zero value of shear stress increase [NO] etc above basal. Tau = 0 checked, changes nothing.
fprintf('current value of tau is %f\n',tau)
rstar=tanh(3.14*tau/lambdacap);                                                         %  GPCR activation as a function of tau, plugs into PIP2 eqns etc in code
W=1.4*(tau/10+sqrt(16*2.86^2+tau^2/100)-4*2.86)^2/(tau/10+sqrt(16*2.86^2+tau^2/100));   % new value of W (the shear-induced calcium current term in calcium eqns) for the new tau above

% calculate the NO production etc etc at the 'new' non-zero shear stress
[T1,X1] = ode15s(@eNOSBiophysJpaperequations_jglfinal_edit,[0 tmax],initValues,options); 

%% tau == 16 dyn/cm^2
tau=16;   % non-zero value of tau, we now check how much a non-zero value of shear stress increase [NO] etc above basal. Tau = 0 checked, changes nothing.
fprintf('current value of tau is %f\n',tau)
rstar=tanh(3.14*tau/lambdacap);                                                         %  GPCR activation as a function of tau, plugs into PIP2 eqns etc in code
W=1.4*(tau/10+sqrt(16*2.86^2+tau^2/100)-4*2.86)^2/(tau/10+sqrt(16*2.86^2+tau^2/100));   % new value of W (the shear-induced calcium current term in calcium eqns) for the new tau above

% calculate the NO production etc etc at the 'new' non-zero shear stress
[T2,X2] = ode15s(@eNOSBiophysJpaperequations_jglfinal_edit,[0 tmax],initValues,options); 

%% tau == 24 dyn/cm^2
tau=24;   % non-zero value of tau, we now check how much a non-zero value of shear stress increase [NO] etc above basal. Tau = 0 checked, changes nothing.
fprintf('current value of tau is %f\n',tau)
rstar=tanh(3.14*tau/lambdacap);                                                         %  GPCR activation as a function of tau, plugs into PIP2 eqns etc in code
W=1.4*(tau/10+sqrt(16*2.86^2+tau^2/100)-4*2.86)^2/(tau/10+sqrt(16*2.86^2+tau^2/100));   % new value of W (the shear-induced calcium current term in calcium eqns) for the new tau above

% calculate signaling dynamics at the 'new' non-zero shear stress
[T3,X3] = ode15s(@eNOSBiophysJpaperequations_jglfinal_edit,[0 tmax],initValues,options); 



%% Plots
figure(1);
for i = 1:4
    subplot(2,2,i); 
    semilogx(T1,X1(:,i),'b-','linewidth',2); hold on
    semilogx(T2,X2(:,i),'g-','linewidth',2); hold on
    semilogx(T3,X3(:,i),'r-','linewidth',2); hold on
    xlabel('Time [s]'); ylabel('Concentration [uM]');
end
suptitle('Concentrations 1-4');
xlabel('Time [s]'); ylabel('Concentration [uM]');

figure(2); 
for i = 1:4
    subplot(2,2,i);
    semilogx(T1,X1(:,i+4),'b-','linewidth',2); hold on
    semilogx(T2,X2(:,i+4),'g-','linewidth',2); hold on
    semilogx(T3,X3(:,i+4),'r-','linewidth',2); hold on
end
suptitle('Concentrations 5-8');
xlabel('Time [s]'); ylabel('Concentration [uM]');

figure(3); 
for i = 1:4
    subplot(2,2,i);
    semilogx(T1,X1(:,i+8),'b-','linewidth',2); hold on
    semilogx(T2,X2(:,i+8),'g-','linewidth',2); hold on
    semilogx(T3,X3(:,i+8),'r-','linewidth',2); hold on
end
suptitle('Concentrations 9-12');
xlabel('Time [s]'); ylabel('Concentration [uM]');

figure(4); 
for i = 1:4
    subplot(2,2,i);
    semilogx(T1,X1(:,i+12),'b-','linewidth',2); hold on
    semilogx(T2,X2(:,i+12),'g-','linewidth',2); hold on
    semilogx(T3,X3(:,i+12),'r-','linewidth',2); hold on
end
suptitle('Concentrations 13-16');
xlabel('Time [s]'); ylabel('Concentration [uM]');
