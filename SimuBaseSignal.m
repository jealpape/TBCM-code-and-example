%close all; clear all; clc;

% addpath('/home/gabriel/Dropbox/MyFunc/VPLab_Code/');
%addpath('PhonationModelsCode')
function [dX,Lg,X,a_g,Ut,Pgo,Pout,Psub,Pcol,tm] = SimuBaseSignal(ActVec,Pl,Vow,gender)

% Simulation information
fs = 44100; % [Hz] sampling frequency
Ns_Larynx = 5;
fs_Larynx = fs/Ns_Larynx;
T_pre = 0.2; N_pre = ceil(T_pre*fs);
T_sim = 0.8; N_Sim = ceil(T_sim*fs);
N_tot = N_pre + N_Sim; t = (1:N_tot)/fs;
N_tran = ceil(T_pre*fs);


% Subglottal pressure
PL = Pl;
% Vowel definition
Vowel = Vow;
  
% Filter design
[b1,a1] = butter(4,50/(44.1e3/2),'high');
[b2,a2] = butter(6,5.5e3/(44.1e3/2));
bf = conv(b1,b2); af = conv(a1,a2);

%% Muscle activation levels  ActVec = [a_LCA, a_IA, a_PCA, a_CT, a_TA]
% ActVec = [a_LCA; a_IA; a_PCA; a_CT; a_TA];
act_cerrada = [0.9;0.9;0;0.5;0.7];
  
%% Larynx object with muscle control
LarynxControl = MuscleControlModel;
LarynxControl.setSimulationParameter(fs_Larynx);
LarynxControl.SelectModelParameters('Alzamendi2020'); % ('Palaparthi2019')  ('Titze2006')  ('Alzamendi2020')

%% Vocal fold model definition and configuration
VFObj = TriangularBodyCoverModel(gender); % For Triangular body-cover vocal fold model
VFObj.setSimulationParameter(fs)
VFObj.setDrivingForceSolver('new');
VFObj.setNonLinMode('on')

%% Vocal tract definition and configuration
VTObj = VocalTractModel(gender);
VTObj.SetSolver('WRA')
if gender(1:4)=='male'
    VTObj.getMaleVocalTract_Story2008(['/' Vowel '/']); % for male
else
    VTObj.getFemaleVocalTract_Story1998(['/' Vowel '/']) % for female
end
VTObj.setSimulationParameter(fs)

%% Subglottal tract definition and configuration
SGTObj = SubglottalTractModel(gender);
SGTObj.SetSolver('WRA')
if strcmpi(Vowel,'NAI')
SGTObj.getSubglottalTract('no-interaction')
else
SGTObj.getSubglottalTract('Story_smooth')
end
VTObj.setSimulationParameter(fs)

% Turbulent flow parameter and initialization
Re_c = 1200;
u_ns_auxb1 = 0;
u_ns_auxb2 = 0;
u_ns_auxa1 = 0;

%% Initializing simulation data struct
X = zeros(N_tot,9);

Pout = zeros(N_tot,1);
Psub = zeros(N_tot,1);
Psup = zeros(N_tot,1);
Pcol = zeros(N_tot,1);

Ut = zeros(N_tot,1);

PGO = zeros(N_tot,1);
a_g = zeros(N_tot,1);
GlotticAngle = nan(N_tot,1);

dX = zeros(N_tot,2); Lg = zeros(N_tot,1);



a_Act = ActVec; % [-] Muscle activation of the five intrinsic muscles

VFObj.setMuscleActivity(a_Act(5),a_Act(4),a_Act(1)-a_Act(3));
LarynxControl.CalcBodyCoverParameters(VFObj);

% Object initialization
VFObj.InitModel
LarynxControl.InitModel
VTObj.InitModel
r_e = 1.0;
SGTObj.InitModel
r_s = 1.0;
  
a_Act = ActVec.*heaviside(t-T_pre)+act_cerrada.*(heaviside(t-0.00)-heaviside(t-T_pre));
% Flow struct information

constFlow = [];
constFlow.PGO = 0; %LarynxControl.aPGO;
constFlow.Ae = VTObj.AreaFunction(1);
constFlow.As = SGTObj.AreaFunction(1);
constFlow.mu = 18.36922e-6;  % Air Viscosity [Pa s]
constFlow.rho = 1.146; % [kg m^-3]  air density
constFlow.c = 350; % speed of sound
constFlow.solver = 'LUCERO2015';
constFlow.L = VFObj.Lg;
constFlow.T = VFObj.Tg;
  
for cont_sim = 1:N_tot
  %% Subglottal pressure transition
  PL_sim = PL * (sin(pi/2*cont_sim/N_tran)*(heaviside(cont_sim)-heaviside(cont_sim-N_tran)) + ...
                 heaviside(cont_sim-N_tran) );

  if rem(cont_sim,Ns_Larynx)==0
    LarynxControl.SimulatePosture(a_Act(:,cont_sim));
    VFObj.setMuscleActivity(a_Act(5,cont_sim),a_Act(4,cont_sim),a_Act(1,cont_sim)-a_Act(3,cont_sim));
    LarynxControl.CalcBodyCoverParameters(VFObj);
  end

  % Simulate Vocal fold dynamics
  Ps_plus = SGTObj.xData(1);
  Pe_minus = VTObj.xData(1);
  VFObj.Simulate(SGTObj.Pressure(1),VTObj.Pressure(1),VTObj.AreaFunction(1))
  a_m = VFObj.ag;

  % Compute glottal volume velocity
  constFlow.PGO = LarynxControl.aPGO;
  constFlow.L = LarynxControl.Lg0 + abs(LarynxControl.Psi_PGO);
  constFlow.T = VFObj.Tg;
  Ug = vfsolver.solveFlow(Ps_plus,Pe_minus,a_m,constFlow);

  % Compute aspiration (noisy) glottal volume velocity (Ref??)
  Re_n = constFlow.rho*Ug/(constFlow.L*constFlow.mu);
  noise_std=1e-12*max(0, Re_n^2-Re_c^2);
  rand_aux = randn;
  u_ns_aux = 1.68*u_ns_auxb1-0.72*u_ns_auxb2+0.05*(rand_aux-0.9*u_ns_auxa1);
  u_ns_auxa1 = rand_aux;
  u_ns_auxb2 = u_ns_auxb1;
  u_ns_auxb1 = u_ns_aux;
  Un = (Re_n - Re_c > 0)*noise_std*u_ns_aux; % Nf(cont_sim)*(Re_n^2 - Re_c^2)*1e-12*(Re_n>Re_c);

  % Total glottal volume velocity
  U = Ug + Un;

  % Simulate acoustic wave propagation throughout subglottal and
  SGTObj.Simulate(U,'PL',PL_sim,r_s)
  VTObj.Simulate(U,r_e) 
  
  
  % Save Quantites 
    X(cont_sim,:) = VFObj.xData';
    PGO(cont_sim) = LarynxControl.aPGO;
    a_g(cont_sim) = a_m + LarynxControl.aPGO;
    Ut(cont_sim) = U;
    
    dX(cont_sim,:) = [max([0, VFObj.xi_02]),max([0, VFObj.xi_01])]; % dxu,dxl
    Lg(cont_sim,:) = VFObj.Lg; 
    
    GlotticAngle(cont_sim) = LarynxControl.angleGlottic;
    Pout(cont_sim) = VTObj.xData(end);
    Psub(cont_sim) = SGTObj.Pressure(1);
    Psup(cont_sim) = VTObj.Pressure(1);
    if VFObj.calcContactArea>0
        [Fu_col,Fl_col] = VFObj.CollisionForces;
        Pcol(cont_sim) = (Fu_col+Fl_col)/VFObj.calcContactArea;
    end

end

tm = round(fs*0.1);


X = 1e3*X(end-tm:end,1:6); % mm

a_g = 1e6*a_g(end-tm:end); %mm 
Ut = 1e6*Ut(end-tm:end); % mL
Pgo = 1e6*PGO(end-tm:end); % mm²

dX = 1e3*dX(end-tm:end,:); % mm
Lg = 1e3*Lg(end-tm:end); % mm


Pout = Pout(end-tm:end);
Psub = Psub(end-tm:end);
Pcol = Pcol(end-tm:end);
    
%%Computing features
% [H1H2,HRF,MFDR,ac_flow,~,SQ,OQ]=get_flow_measures_Journal_Z_Scores(1e3*DataSim.Ut(end-n_rest:end)',fs);
% F0 = measure.f0(DataSim.Ut(end-n_rest:end),fs);
% Pcol = measure.max(DataSim.ContactP(end-n_rest:end),fs);
% Psub = mean(DataSim.Psub(end-n_rest:end));
% SPL = measures_getSPL(DataSim.Pout(end-n_rest:end),fs);

end
