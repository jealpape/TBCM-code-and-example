close all; clear all; clc;

% addpath('/home/gabriel/Dropbox/MyFunc/VPLab_Code/');
% addpath('../zz_PhonationModel_MatlabClasses/')

% Simulation information
fs = 44100; % [Hz] sampling frequency
Ns_Larynx = 5;
fs_Larynx = fs/Ns_Larynx;
T_pre = 0.4; N_pre = ceil(T_pre*fs);
T_sim = 0.3; N_Sim = ceil(T_sim*fs);
N_tot = N_pre + N_Sim; t = (1:N_tot)/fs;
N_tran = ceil(0.1*T_pre*fs);

% Subglottal pressure
% PL = (0.8+0.6*a_CT+0.4*a_TA)*1e3; 
PL = 800;

% Vocal tract cases to investigate
VowelCases = {'A','i'};%,'AE','I','e','u','o','U','E','V','O'};
L_VowelCases = length(VowelCases);

for cont_VowelCase = 1:L_VowelCases
  Vowel = VowelCases{cont_VowelCase};
  if any(strcmp(Vowel(1),{'a' 'e' 'i' 'o' 'u'}))
    flag = '0';
  else
    flag = '1';
  end
  
  % Filter design
  [b1,a1] = butter(4,50/(44.1e3/2),'high');
  [b2,a2] = butter(6,5.5e3/(44.1e3/2));
  bf = conv(b1,b2); af = conv(a1,a2);

  % Muscle activation levels  ActVec = [a_LCA, a_IA, a_PCA, a_CT, a_TA]
  ActVec = [0.50, 0.50, 0.0, 0.1, 0.25];

  % Larynx object with muscle control
  LarynxControl = MuscleControlModel;
  LarynxControl.setSimulationParameter(fs_Larynx);
  LarynxControl.SelectModelParameters('Alzamendi2020'); % ('Palaparthi2019')  ('Titze2006')  ('Alzamendi2020')

  % Vocal fold model definition and configuration
%   VFObj = BodyCoverModel('male'); % For Body-cover vocal fold model
  VFObj = TriangularBodyCoverModel('male'); % For Triangular body-cover vocal fold model
  VFObj.setSimulationParameter(fs)
  VFObj.setDrivingForceSolver('new');

  % Vocal tract definition and configuration
  VTObj = VocalTractModel('male');
  VTObj.SetSolver('WRA')
  if strcmpi(Vowel,'NAI')
    VTObj.getSimpleVocalTract('no-interaction')
  else
    VTObj.getMaleVocalTract_Story2008(['/' Vowel '/']) % For male voice
%     VTObj.getFemaleVocalTract_Story1998(['/' Vowel '/']) % For female voice
  end
  VTObj.setSimulationParameter(fs)

  % Subglottal tract definition and configuration
  SGTObj = SubglottalTractModel('male');
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

  % Initializing simulation data struct
  DataSim.X_masses = zeros(N_tot,9);
  DataSim.Pout = zeros(N_tot,1);
  DataSim.Psub = zeros(N_tot,1);
  DataSim.Psup = zeros(N_tot,1);
  DataSim.audiosig = zeros(N_tot,1);
  DataSim.Ug = zeros(N_tot,1);
  DataSim.Ut = zeros(N_tot,1);
  DataSim.Uout = zeros(N_tot,1);
  DataSim.a_m = zeros(N_tot,1);
  DataSim.a_g = zeros(N_tot,1);
  DataSim.StrainVF = nan(N_tot,1);
  DataSim.GlotticAngle = nan(N_tot,1);
  DataSim.ContactP = zeros(N_tot,1);
  DataSim.a_Act = zeros(5,1);
  DataSim.fs = fs;
  DataSim.PL = PL;
  DataSim.Ns_Larynx = Ns_Larynx;

    a_Act = ActVec; % [-] Muscle activation of the five intrinsic muscles
    DataSim.a_Act = a_Act;

    % Object initialization
    LarynxControl.CalcBodyCoverParameters(VFObj,a_Act(5));
    VFObj.InitModel
    LarynxControl.InitModel
    VTObj.InitModel
    r_e = 1.0;
    SGTObj.InitModel
    r_s = 1.0;
  
    % Flow struct information
    constFlow = [];
    constFlow.PGO = LarynxControl.aPGO;
    constFlow.Ae = VTObj.AreaFunction(1);
    constFlow.As = SGTObj.AreaFunction(1);
    constFlow.mu = 18.36922e-6;  % Air Viscosity [Pa s]
    constFlow.rho = 1.146; % [kg m^-3]  air density
    constFlow.c = 350; % speed of sound
    constFlow.solver = 'LUCERO2015';
    constFlow.L = VFObj.Lg;
    constFlow.T = VFObj.Tg;
  
    for cont_sim = 1:N_tot
      % Subglottal pressure transition
      PL_sim = PL * (sin(pi/2*cont_sim/N_tran)*(heaviside(cont_sim)-heaviside(cont_sim-N_tran)) + ...
                     heaviside(cont_sim-N_tran) );

      if rem(cont_sim,Ns_Larynx)==0
        LarynxControl.SimulatePosture(a_Act);
        LarynxControl.CalcBodyCoverParameters(VFObj,a_Act(5));
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
    
      % Compute aspiration (noisy) glottal volume velocity
      Re_n = constFlow.rho*Ug/(constFlow.L*constFlow.mu);
      noise_std=1e-12*max(0, Re_n^2-Re_c^2);
      rand_aux = randn;
%       u_ns_aux = 1.4*u_ns_auxb1-0.49*u_ns_auxb2+0.05*(rand_aux-0.8*u_ns_auxa1);
      u_ns_aux = 1.68*u_ns_auxb1-0.72*u_ns_auxb2+0.05*(rand_aux-0.9*u_ns_auxa1);
      u_ns_auxa1 = rand_aux;
      u_ns_auxb2 = u_ns_auxb1;
      u_ns_auxb1 = u_ns_aux;
      Un = (Re_n - Re_c > 0)*noise_std*u_ns_aux; % Nf(cont_sim)*(Re_n^2 - Re_c^2)*1e-12*(Re_n>Re_c);
    
      % Total glottal volume velocity
      Ut = Ug + Un;
    
      % Simulate acoustic wave propagation throughout subglottal and
      % supraglottal tracts.
      SGTObj.Simulate(Ut,'PL',PL_sim,r_s)
      VTObj.Simulate(Ut,r_e) 
    
        DataSim.X_masses(cont_sim,:) = VFObj.xData';
        DataSim.a_m(cont_sim) = a_m;
        DataSim.a_g(cont_sim) = a_m + LarynxControl.aPGO;
        DataSim.Ug(cont_sim) = Ug;
        DataSim.Ut(cont_sim) = Ut;
        DataSim.StrainVF(cont_sim) = LarynxControl.getStrainVF;
        DataSim.GlotticAngle(cont_sim) = LarynxControl.angleGlottic;
        DataSim.Pout(cont_sim) = VTObj.xData(end);
        DataSim.Uout(cont_sim) = VTObj.Airflow(end);
        DataSim.Psub(cont_sim) = SGTObj.Pressure(1);
        DataSim.Psup(cont_sim) = VTObj.Pressure(1);
        if VFObj.calcContactArea>0
            [Fu_col,Fl_col] = VFObj.CollisionForces;
            DataSim.ContactP(cont_sim) = (Fu_col+Fl_col)/VFObj.calcContactArea;
        end
        
    end
    
    filename = ['Sim_SustainedPhonation_' Vowel flag '.wav'];
    audiosig = filtfilt(bf,af,DataSim.Pout(end-N_Sim:end));
    audiosig = 0.9*audiosig/max(abs(audiosig));
    audiosig = audiosig.*tukeywin(length(audiosig),0.02);
    DataSim.audiosig = audiosig;
    audiowrite(filename,audiosig,fs)

  
  figure
  subplot(4,1,1)
  plot(t,DataSim.X_masses(:,1:3)*1e3)
  ylabel('Disp. in mm')
  xlabel('Time in s')
  legend('x_u','x_l','x_b')
  title('Mass displacement')
  
  subplot(4,1,2)
  plot(t,DataSim.a_g*1e4)
  ylabel('a_{g} in cm^2')
  xlabel('Time in s')
  title('Glottal area')
  
  subplot(4,1,3)
  plot(t,DataSim.Ut*1e6)
  ylabel('U_{g} in ml/s')
  xlabel('Time in s')
  title('Glottal volume velocity')
  
  subplot(4,1,4)
  plot(t,DataSim.Pout)
  ylabel('P_{o} in Pa')
  xlabel('Time in s')
  title('Radiated voice pressure')

  % Erasing result structure and objects
  DataSim = [];
  LarynxControl = [];
  VFObj = [];
  SGTObj = [];
  VTObj = [];

end