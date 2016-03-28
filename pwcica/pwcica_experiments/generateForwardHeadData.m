%% (OPTIONAL, need SIFT)  Single-trial, Time-Varying Simulation (Mullen, 2011 unpublished)
%  System of stocastically-forced, damped coupled oscillators with
%  time-varying coupling coefficients
%  This simulation creates a simualted "seizure" with time-varying coupling
%  between clusters of sources which switches between 4 different stages.
%PWC-ICA Paper; 
% Experiment 4 (Forward model);
% NOTE: Requires EEGLAB
% Instructions:
% Either load pwcica_exp4.mat, or if you have the SIFT toolbox you may run
% the code in the first section to generate a randomized version of the
% data.

%IMPORTANT: Set head model paths:
HeadModelPath = [pwd,'\Colin27_Biosemi_1010_standard.mat'];
forwardParametersPath = [pwd,'forwardParameters.mat'];

SamplingRate = 100; % Hz

Nl = 5*60*SamplingRate;  % length of each epoch (5 minutes)
Nr = 1;                  % number of trials (realisations)
ndisc = 1000;            % number of samples to discard from VAR model (startup transients)
ModelOrder = 6;          % VAR model order

% Set the fundamental frequency (f) and damping time (tau) of the 
% oscillators for each cluster of sources
% NOTE: A small tau generates a "noisier" signal

% Cluster 1
f1=20; 
tau1 = 20;

% Cluster 2
f2=20;
tau2 = 3;

% Cluster 3
f3=10;
tau3 = 7;

% Cluster 4
f4 = 10;
tau4 = 6;


% set the approximate durations of each stage
S1_width = 5*SamplingRate;
S2_width = 5*SamplingRate;
S3_width = 5*SamplingRate;
S4_width = 5*SamplingRate;

S0 = 60*SamplingRate; % start of seizure (1 minute)

S1_center = S0+S1_width/2;                      % center of stage 1
S2_center = S1_center+S1_width;                 % center of stage 2
S3_center = S2_center;                          % center of stage 3
S4_center = S3_center+S3_width;                 % center of stage 4

Offset = 0;

% write the system of equations for coupled damped oscillators
expr = { ...
      sprintf('x1(t) = {2*exp(-1/(%f+normpdfg(t,%f,%f,%f,100)))*cos(2*pi*20.000000/100.000000)}*x1(t-1) + {-exp(-2/(%f+normpdfg(t,%f,%f,%f,100)))}*x1(t-2) + e1(t)',tau1,S1_width/2,8,Offset+S1_center, tau1,S1_width,8,Offset+S1_center) ...                    % Ictal driver
      sprintf('x2(t) = %s + {normpdfg(t,%f,%f,%f,0.01)}*x3(t-2) + {normpdfg(t,%f,%f,%f,0.01)}*x4(t-2)  + {normpdfg(t,%f,%f,%f,1.3)}*x1(t-6)    + e2(t)',sim_dampedOscillator(f2,tau2,SamplingRate,2),  S1_width/2,8,Offset+S1_center,   S1_width/2,8,Offset+S1_center,  S1_width/2,8,Offset+S1_center) ...                     % CLUSTER 2
      sprintf('x3(t) = %s +  {normpdfg(t,%f,%f,%f,0.7)}*x2(t-2) + {normpdfg(t,%f,%f,%f,0.1)}*x4(t-2)  + {normpdfg(t,%f,%f,%f,0.3)}*x5(t-3)    + e3(t)',sim_dampedOscillator(f2-1,tau2,SamplingRate,3),  S1_width/2,8,Offset+S1_center,   S1_width/2,8,Offset+S1_center,   (S2_width+S3_width)/2,10,Offset+S3_center) ...       % CLUSTER 2
      sprintf('x4(t) = %s +  {normpdfg(t,%f,%f,%f,0.7)}*x2(t-2) + {normpdfg(t,%f,%f,%f,0.1)}*x3(t-2)                                          + e4(t)',sim_dampedOscillator(f2+1,tau2,SamplingRate,4),  S1_width/2,8,Offset+S1_center,  S1_width/2,8,Offset+S1_center) ... % CLUSTER 2
      sprintf('x5(t) = %s +  {normpdfg(t,%f,%f,%f,0.001)}*x6(t-2) + {normpdfg(t,%f,%f,%f,0.5)}*x3(t-3)        + e5(t)'  ,sim_dampedOscillator(f3,tau3,SamplingRate,5), S3_width/2,8,Offset+S3_center , S3_width/2,8,Offset+S3_center) ...  % CLUSTER 3
      sprintf('x6(t) = %s +  {normpdfg(t,%f,%f,%f,0.001)}*x5(t-2)                                             + e6(t)'  ,sim_dampedOscillator(f3,tau3,SamplingRate,6), S3_width/2,8,Offset+S3_center) ...                                  % CLUSTER 3
      sprintf('x7(t) = %s + {normpdfg(t,%f,%f,%f,1.3)}*x4(t-6) + {normpdfg(t,%f,%f,%f,0.01)}*x9(t-2) + {normpdfg(t,%f,%f,%f,0.01)}*x8(t-2) + {normpdfg(t,%f,%f,%f,0.01)}*x10(t-2)    + e7(t)' ,sim_dampedOscillator(f4,tau4,SamplingRate,7),     S4_width/2,8,Offset+S4_center,     S4_width/2,8,Offset+S4_center,     S4_width/2,8,Offset+S4_center,     S4_width/2,8,Offset+S4_center) ...  % CLUSTER 4
      sprintf('x8(t) = %s + {normpdfg(t,%f,%f,%f,0.6)}*x7(t-2)         + e8(t)' ,sim_dampedOscillator(f4,tau4,SamplingRate,8),     S4_width/2,8,Offset+S4_center) ... % CLUSTER 4
      sprintf('x9(t) = %s + {normpdfg(t,%f,%f,%f,0.6)}*x7(t-2)         + e9(t)' ,sim_dampedOscillator(f4+1,tau4,SamplingRate,9),   S4_width/2,8,Offset+S4_center) ... % CLUSTER 4
      sprintf('x10(t)= %s + {normpdfg(t,%f,%f,%f,0.6)}*x7(t-2)         + e10(t)',sim_dampedOscillator(f4-1,tau4,SamplingRate,10),  S4_width/2,8,Offset+S4_center) ... % CLUSTER 4
      sprintf('x11(t)= 1.3*x11(t-1) + -0.8*x11(t-2)                             + e11(t)') ...    % Decoupled
      sprintf('x12(t)= 1.2*x12(t-1) + -0.4*x12(t-2)                             + e12(t)') ...    % Decoupled
      sprintf('x13(t)= 0.8*x13(t-1) + -0.4*x13(t-2) + -0.1*x13(t-4)             + e13(t)') ...    % Decoupled
      };
  
% create prototype VAR structure
Aproto = sim_genVARModelFromEq(expr,ModelOrder);

% STEP 2: Simulate the VAR process

[A] = sim_genTVARcoeffs(Aproto,ModelOrder,Nl,'NumSamplesToDiscard',ndisc,'Verbose',true);

% STEP 3: generate data from the VAR model

% Specify the noise covariance matrix. 
% Sigma is the noise variance.
sigma = 1;
M = size(A{1},1);
C = sigma*eye(M);             

% % hyperbolic secant noise
% data = permute(tvarsim(zeros(1,M),A,C,[Nl Nr],ndisc,2/pi,0,'hsec'),[2 1 3]);

% laplacian noise (generalized gaussian)
% data = permute(tvarsim(zeros(1,M),A,C,[Nl Nr],ndisc,1,1,'gengauss'),[2 1 3]);

% gaussian noise
data = permute(tvarsim(zeros(1,M),A,C,[Nl Nr],ndisc,2,1,'gengauss'),[2 1 3]);


% STEP 4: Create EEGLAB dataset and source potentials
% 
EEG = eeg_emptyset;
EEG.data = data;
[EEG.icaweights EEG.icasphere EEG.icawinv] = deal(eye(M));
EEG.icaact = [];
EEG.srate = SamplingRate;
EEG.times = ((0:(Nl-1))/SamplingRate)*1000;   % ms
EEG.pnts = Nl;
EEG.trials = Nr;
EEG.xmin = EEG.times(1);
EEG.xmax = EEG.times(end)/1000;  % sec
EEG.nbchan = M;
EEG.setname = 'VAR Simulation';
EEG.condition = 'VAR Simulation';
EEG.srcpot = data;

%
EEG = eeg_checkset(EEG);
pop_eegplot(EEG);

ALLEEG = EEG;
CURRENTSET = length(ALLEEG);
eeglab redraw;

% Do forward model:
%[hlp_getSiftRoot filesep 'resources' filesep 'headmodels' filesep 'standard-Colin27-385ch.mat'];
hmObj = hlp_validateHeadModelObject(HeadModelPath);

roiNames = hmObj.atlas.label;
% roiIndex = randperm(length(roiNames),nsrc);
% roiNames = roiNames(roiIndex);

% roiIndex = [16    23    11    15    44     8    88    35    89     6    65    67    45];
roiIndex = [122    16   115    87     9   135   145     2   146    41    49    26   144];
roiNames = roiNames(roiIndex);

% project data through forward model
% ----------------------------------------------------------------------------------------------
[EEGmerge, ~, fwdModel] = sim_eegdata('srcdyn',{'arg_selection' 'Precomputed' 'srcdata' data}, ... %'Precomputed'
                                        'fwdproj', ...
                                            {'hmObj' hmObj ...
                                            'sourceAtlasLabels' {} ...
                                            'LFM' []   ...
                                            'channels' {}  ...
                                            'sourceShape' {'arg_selection' 'dipole' 'roiAtlasLabels' {} 'roiOrdered' roiNames 'nearestNeighbor' true}, ...
                                            'addNoise', {'SignalToNoise' 4}, 'verb' 2, }, ...
                                        'makedipfit',true,'vismodel',false,'verb',2);

% cleanup
% hmObj.saveToFile(segments.HeadModelFileName);
% delete(hmObj.surfaces);
% delete(hmObj.leadFieldFile);
        
EEGmerge.srate   = SamplingRate;
EEGmerge.pnts    = size(EEGmerge.data,2);
EEGmerge.times   = 0:(EEGmerge.pnts-1);%linspace(0,(EEGmerge.pnts)/SamplingRate,EEGmerge.pnts);
EEGmerge.xmin    = EEGmerge.times(1)/SamplingRate;
EEGmerge.xmax    = EEGmerge.times(end)/SamplingRate;
EEGmerge.trials = 1;
% store the lead-field matrix
for k=1:length(EEG)
    EEGmerge.etc.LFM{k} = fwdModel(k).LFM(:,fwdModel(k).centroids_LFM_idx);
    % preserve original channel locations
    EEGmerge.etc.chanlocs_ori{k} = EEG(k).chanlocs;
end
EEGmerge.chaninfo.nosedir = '+X';
EEGmerge.chaninfo.plotrad = 0.5;
EEGmerge.etc.roiNames = roiNames;
EEGmerge.etc.roiIndex = roiIndex;
EEGmerge = rmfield(EEGmerge,'srcpot_all');
% assign the setname based on head model name
%[~,pm,~] = fileparts(segments(1).HeadModelFileName);
EEGmerge.setname = ['VAR Sim Data '];% pm];