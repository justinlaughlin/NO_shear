%%
function dy = f(t,y)
dy = zeros(16,1);    % a column vector

% global variables
global tau mu1 rr aktt pkct kcatp camtot bufftot pip2t phiinf kleak delta Kcoup SGC0 kpip2 kmpip2 kmakt kcam0 kmcam keakt kmeakt kpt5cam qmax
global kcce cs0 cex kres krel kout vr VpPKC kcavdec Gt ka kd alp api3k rstar W rcgmp Vdg kmdg NOdec RNO eta k2 k3 k5 k2p theta N xi kmpkc kdcam
global  kca4cam kon a0 a1 g0 g1 b1

G               = y(1);
IP3             = y(2);
PIP2            = y(3);
PIP3            = y(4);
CAc             = y(5);
CAs             = y(6);
CAb             = y(7);
AKTa            = y(8);
PKCa            = y(9);
Ca4CaM          = y(10);
eNOS_CaM        = y(11);
eNOS_CaMa       = y(12);
eNOScav0        = y(13);
NO              = y(14);
cGMP            = y(15);
eNOScav         = y(16);

% Hill coefficient in Ca4CaM equation
Beta    = 2.7;

% Intermediate Equations (JGL)
Psi     = 1.00 - xi*cGMP;                          % (10). called by (9)
Gd = Gt-G;

rf = alp*Kcoup*phiinf/(phiinf+kcatp)*G; 
kp_pip2 = kpip2/(1 + api3k)* (1+api3k*exp(-eta*t)*tanh(pi*tau/delta)) + k2p;   % (12) activation of PI3K increases rate of PIP2 phosphorylation to PIP3 

qrel    = krel*(IP3/(k2+IP3))^3*CAs;                           	 % (3) calcium release from internal stores
qres    = kres*(CAc/(k3+CAc))^2;                               	 % (4) calcium re-sequestration into internal stores
qout    = kout*(CAc/(k5+CAc));                         	         % (5) calcium efflux via sodium-calcium exchanger             
qmsic   = qmax/(1+N*exp(-W));                       	 	 % (7) 
qcce    = kcce*Psi*(cs0 - CAs)*(cex - CAc);             	 % (9) 
qin     = qmsic + qcce;
qon     = kon*CAc*(bufftot-CAb);
qoff    = koff*CAb;


PIP3p   = PIP3/pip2t*100;				       	% PIP3 percentage 
kpakt   = .1*kmakt*(PIP3p-.31)/(3.1-.31);
kppkc   = .1*kmpkc*(PIP3p-.31)/(3.1-.31);
PKC     = pkct - PKCa;
AKT     = aktt - AKTa;

% ODEs
dy(1)= ka*Gd*rstar - kd*G;
dy(2) = rf*PIP2 - mu1*IP3;
dy(3) = -(rf+rr)*PIP2 - rr*IP3 + rr*pip2t-(PIP2*kp_pip2 - kmpip2*PIP3);
dy(4) = PIP2*kp_pip2 - kmpip2*PIP3; 
dy(5) = qrel - qres - qout + qin + kleak*(CAs)^2 - qon + qoff;
dy(6) = -vr*(qrel - qres +kleak*(CAs)^2);
dy(7) = qon - qoff;
dy(8) = kpakt*AKT - kmakt*AKTa; 
dy(9) = kppkc*PKC - kmpkc*PKCa;
dy(10) = kca4cam*(camtot*theta*CAc^Beta/(kdcam+CAc^Beta) - Ca4CaM);   
dy(11) = kcam0*Ca4CaM/(kpt5cam+Ca4CaM)*eNOScav - kmcam*eNOS_CaM - (keakt*eNOS_CaM*AKTa/aktt - kmeakt*eNOS_CaMa); 
dy(12) = keakt*eNOS_CaM*AKTa/aktt - kmeakt*eNOS_CaMa; 
dy(13) = VpPKC*PKCa/pkct*(eNOScav) - kcavdec*eNOScav0; 
dy(14) = RNO*(eNOS_CaM+9*eNOS_CaMa) - NOdec*NO - 0.022*SGC0*(NO^2 + NO*b1)/(a0+a1*NO+NO^2);
dy(15) = rcgmp*(g0+g1*NO+(NO)^2)/(a0+a1*NO+(NO)^2) - (cGMP)^2*Vdg/(kmdg+cGMP);
dy(16) = - (VpPKC*PKCa/pkct*(eNOScav) - kcavdec*eNOScav0) - (kcam0*Ca4CaM/(kpt5cam+Ca4CaM)*eNOScav - kmcam*eNOS_CaM);     

end
