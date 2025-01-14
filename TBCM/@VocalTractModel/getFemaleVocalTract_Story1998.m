%%
% getFemaleVocalTract_Story1998: Function to get the female vocal tract
% configurations reported in [1], and the relevant simulation parameters,
% in order to configure the vocal tract model accordingly.
%
% Structure: getFemaleVocalTract_Story1998(VTobj,VTConf)
%            Data = getFemaleVocalTract_Story1998(VTobj,VTConf)
%            [AreaFun,Fsamp] = getFemaleVocalTract_Story1998(VTobj,VTConf)
%            [AreaFun,Fsamp,Deltaz] = getFemaleVocalTract_Story1998(VTobj,VTConf)
%
% where
%
% VTobj: is an object from VocalTractModel (handle) class,
% VTConf: string for the selected idealized vocal tract configuration. 
%         Available options = {'/i/','/I/','/E/','/AE/','/V/','/A/','/O/',
%                              '/o/', '/U/','/u/','/3/' }.
% Data: Struct gathering the relevan parameter for the selected vocal tract
%       configuration.
% AreaFun: Area function for the selected vocal tract configuration.
% Fsamp: Sampling frequency.
% Deltaz: length if the vocal tract sections.
%
% References:
% [1] B. H. Story, I. R. Titze, and E. A. Hoffman, “Vocal tract area functions 
% for an adult female speaker based on volumetric imaging,” J. Acoust. Soc. Am., 
% vol. 104, no. 1, p. 471, Jun. 1998.
%
% Coded by Gabriel Alzamendi, October 2020.
% Based on original code by Matías Zañartu and Gabriel Galindo

function varargout = getFemaleVocalTract_Story1998(VTobj,varargin)
  VOCALTRACTOPTIONS = {'/i/','/I/','/E/','/AE/','/V/','/A/','/O/', ...
                       '/o/','/U/','/u/','/i_ct/','/A_ct/' };
  SampFreq = 44.1e3; % [Hz] Sampling frequency for simulating the idealized
                     %  vocal tract using WRA.
  LenghtSec = 3.96825e-3; % [m] Length of vocal tract sections
  if (nargin == 1)
    tract_type = '/A/';
  elseif (nargin == 2)&&any(strcmp(varargin{1},VOCALTRACTOPTIONS))
    tract_type = char(varargin{1});
  else
    SolverOpts = '';
    for OptsAux = VOCALTRACTOPTIONS
      SolverOpts = strcat(SolverOpts, ' ''', char(OptsAux), ''', ');
    end
    SolverOpts = SolverOpts(1:end-1);
    error('Vocal tract configuration not available! The acceptable options are:%s. \n',SolverOpts);  
  end
  % Available idealized vocal tract configurations
  switch tract_type
    case '/i/'  % Close front unrounded vowel: IPA number 301
      section = 10^-4*[
         0.36 0.50 0.82 1.00 1.45 3.90 3.33 2.84 2.75 3.84 ...
         4.30 4.38 4.16 3.22 2.50 2.25 2.19 2.27 1.78 1.25 ...
         0.74 0.28 0.18 0.14 0.19 0.28 0.48 0.75 0.51 0.71 ...
         1.00 0.74 1.24 1.51 1.65 1.59 ];
         
    case '/i_ct/'  % Close front unrounded vowel: IPA number 301
      section = 10^-4*[
         0.50 0.36 0.38 0.49 1.45 2.03 1.87 1.94 2.13 3.02 ...
         3.90 3.73 3.05 2.27 1.77 1.73 1.60 1.21 0.68 0.46 ...
         0.21 0.03 0.04 0.05 0.07 0.14 0.27 0.32 0.41 0.53 ...
         0.69 1.14 1.80 2.34 2.01 1.95 ];
        
    case '/I/'  % Near-close near-front unrounded vowel: IPA number 319
      section = 10^-4*[
        0.21 0.14 0.40 0.51 0.82 0.72 0.88 1.01 0.91 1.55 ...
        1.70 1.86 1.36 1.14 1.67 1.93 1.79 1.57 1.46 1.22 ...
        1.28 1.22 1.31 1.23 1.13 1.20 0.98 1.47 1.70 0.97 ...
        0.86 1.20 1.18 1.34 ];
      
    case '/E/'  % Open-mid front unrounded vowel: IPA number 303
      section = 10^-4*[
        0.18 0.15 0.55 0.92 0.49 0.53 0.50 0.46 1.35 1.16 ...
        1.54 1.31 1.03 1.38 1.67 1.45 1.26 1.22 1.29 1.27 ...
        1.28 1.45 1.83 1.46 1.26 1.25 1.77 1.60 1.52 1.80 ...
        2.14 1.59 ];
      
    case '/AE/'  % Near-open front unrounded vowel: IPA number 325
      section = 10^-4*[
        0.26 0.19 0.70 1.23 1.14 1.03 1.10 1.27 1.78 1.94 ...
        1.80 2.16 2.78 3.01 3.00 2.98 2.92 3.27 3.80 3.79 ...
        4.03 3.46 3.06 2.75 3.89 2.91 3.22 3.82 3.90 3.43 ];
      
    case '/V/'  % Open-mid back unrounded vowel: IPA number 314
      section = 10^-4*[
        0.18 0.29 0.51 0.94 0.43 0.37 0.42 0.27 0.51 0.60 ...
        1.00 0.48 0.21 0.28 0.62 0.84 0.92 1.15 1.47 2.00 ...
        2.32 2.71 2.64 2.74 2.26 2.14 1.96 2.02 1.83 1.25 ...
        1.35 1.90 1.70 1.69 ];
      
    case '/A/'  % Open back unrounded vowel: IPA number 305
      section = 10^-4*[
         0.20 0.42 0.46 1.17 1.55 0.99 0.38 0.18 0.97 1.17 ...
         1.24 1.32 1.08 1.18 1.38 2.04 2.06 2.36 2.78 2.86 ...
         3.50 3.54 3.84 4.31 5.16 5.19 4.39 3.20 2.89 2.13 ...
         1.95 1.38 1.29 2.81 ];
      
    case '/A_ct/'  % Open back unrounded vowel: IPA number 305
      section = 10^-4*[
         0.56 0.24 0.19 0.20 0.46 0.80 0.47 0.25 0.22 0.30 ...
         0.42 0.75 1.40 1.01 0.78 1.26 1.54 1.64 1.98 2.40 ...
         3.16 3.65 4.33 4.35 4.58 4.94 4.59 4.39 3.98 3.93 ...
         3.72 3.37 2.81 2.81 3.11 2.94 ];
      
    case '/O/'  % Open-mid back rounded vowel: IPA number 306
      section = 10^-4*[
        0.82 1.11 0.95 0.70 0.44 0.45 0.38 0.36 0.56 1.01 ...
        0.95 0.79 0.68 0.33 0.88 2.06 2.41 2.47 3.10 4.01 ...
        4.43 5.15 5.61 6.07 6.62 7.15 8.05 8.38 8.13 8.10 ...
        8.64 6.02 3.32 2.26 1.94 1.58 ];
      
    case '/o/'  % Close-mid back rounded vowel: IPA number 307
      section = 10^-4*[
        0.17 0.18 0.81 0.65 0.61 0.46 0.44 0.61 0.67 0.97 ...
        1.19 0.66 0.37 0.46 1.33 2.00 1.81 1.77 2.24 2.43 ...
        2.87 4.24 4.63 4.95 4.92 4.76 5.34 5.43 4.75 4.26 ...
        2.97 1.99 1.33 1.45 1.60 1.39 ];
      
    case '/U/'  % Near-close near-back rounded vowel: IPA number 321
      section = 10^-4*[
        0.07 0.10 0.07 0.05 0.18 0.43 0.75 0.47 0.29 0.11 ...
        0.32 0.58 1.11 0.84 0.17 0.39 0.81 1.03 1.02 1.04 ...
        1.42 1.93 2.10 2.21 2.53 2.82 2.80 2.46 2.32 2.46 ...
        2.09 2.57 1.76 1.37 1.09 1.16 ];
      
    case '/u/'  % Close back rounded vowel: IPA number 308
      section = 10^-4*[
        0.22 0.48 0.81 0.90 1.20 3.13 2.71 1.75 2.01 2.20 ...
        2.16 2.56 2.61 2.34 2.05 1.21 0.80 0.74 1.36 1.18 ...
        1.05 0.78 0.78 1.17 1.47 2.00 2.50 2.92 3.10 2.98 ...
        3.40 3.61 2.97 2.74 1.79 0.87 0.44 0.18 ];
      
  end
  section = section(:);
  N_sections = length(section);
  
  % Gender checking
  if ~strcmp(VTobj.sex,'female')
      warning('Sex setting is changed to ''female''!');
  end
  % Setting model parameters
  VTobj.AreaFunction = section;
  VTobj.N_AreaSection = N_sections;
  VTobj.Delta_z = LenghtSec;
  VTobj.sex = 'female';
  VTobj.setSimulationParameter(SampFreq);
  
  % Function Output
  switch nargout
    case 1
      Info.AreaFunction = section;
      Info.N_AreaSection = N_sections;
      Info.Delta_z = LenghtSec;
      Info.sex = 'female';
      Info.fs = SampFreq;
      Info.Ts = 1/SampFreq;
      varargout{1} = Info;
    case 2
      varargout{1} = section;
      varargout{2} = SampFreq;
    case 3
      varargout{1} = section;
      varargout{2} = SampFreq;
      varargout{3} = LenghtSec;
  end
  
end
