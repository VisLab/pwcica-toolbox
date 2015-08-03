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

load('pwcica_exp4.mat');
%% (OPTIONAL, need SIFT)  Single-trial, Time-Varying Simulation (Mullen, 2011 unpublished)
%  System of stocastically-forced, damped coupled oscillators with
%  time-varying coupling coefficients
%  This simulation creates a simualted "seizure" with time-varying coupling
%  between clusters of sources which switches between 4 different stages.

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

%% Run infomax

tic;
[EEGim] = pop_runica(EEGmerge,'extended',1);
timeElapsed.im = toc;
Wim = EEGim.icaweights*EEGim.icasphere;

%% Run FastICA

tic;
[~,Wfast] = fastica(EEGmerge.data,'approach','symm');
timeElapsed.fast = toc;

%% Run PWCICA

tic;
[W1,metas1] = pwcica(EEGmerge,'Level',1,'TimeInvariant',0,'NoPaths',1);
timeElapsed.pwcica1 = toc;
tic;
[W2,metas2] = pwcica(EEGmerge,'Level',2,'TimeInvariant',0,'NoPaths',1);
timeElapsed.pwcica2 = toc;
tic;
[W4,metas4] = pwcica(EEGmerge,'Level',4,'TimeInvariant',0,'NoPaths',1);
timeElapsed.pwcica4 = toc;
tic;
[W8,metas8] = pwcica(EEGmerge,'Level',8,'TimeInvariant',0,'NoPaths',1);
timeElapsed.pwcica8 = toc;

%% Run PWCICA with Haar

tic;
[W1h,metas1h] = pwcica(EEGmerge,'Level',1,'TimeInvariant',1,'NoPaths',1);
% W1h = (real(metas1h.complexW*metas1h.SphC)-4*imag(metas1h.complexW*metas1h.SphC))*metas1h.Sph;
timeElapsed.pwcica1h = toc;
tic;
[W2h,metas2h] = pwcica(EEGmerge,'Level',2,'TimeInvariant',1,'NoPaths',1);
% W2h = (real(metas2h.complexW*metas2h.SphC)-4*imag(metas2h.complexW*metas2h.SphC))*metas2h.Sph;
timeElapsed.pwcica2h = toc;
tic;
[W4h,metas4h] = pwcica(EEGmerge,'Level',4,'TimeInvariant',1,'NoPaths',1);
% W4h = (real(metas4h.complexW*metas4h.SphC)-4*imag(metas4h.complexW*metas4h.SphC))*metas4h.Sph;
timeElapsed.pwcica4h = toc;
tic;
[W8h,metas8h] = pwcica(EEGmerge,'Level',8,'TimeInvariant',1,'NoPaths',1);
timeElapsed.pwcica8h = toc;
tic;
[W64h,metas64h] = pwcica(EEGmerge,'Level',64,'TimeInvariant',1,'NoPaths',10);
timeElapsed.pwcica64h = toc;

%% Run Hilbert transform approach
tic;
[WH,metasH] = pwcica(EEGmerge,'Level',0,'TimeInvariant',1,'NoPaths',1);
% W1h = (real(metas1h.complexW*metas1h.SphC)-4*imag(metas1h.complexW*metas1h.SphC))*metas1h.Sph;
timeElapsed.pwcicaH = toc;

%% Non-greedy
EEG.srcpot = data;
figure; hold on;
p(1) = errorbar([1:13],...
    mean(abs(corr(EEGmerge.srcpot',(EEGmerge.data)')')),...
    std(abs(corr(EEGmerge.srcpot',(EEGmerge.data)')')),...
    std(abs(corr(EEGmerge.srcpot',(EEGmerge.data)')')),'g','LineWidth',2);

p(2) = plot( max(abs(corr(EEG.srcpot',(WH*EEGmerge.data)')')),'c','LineWidth',3);
p(3) = plot( max(abs(corr(EEG.srcpot',(Wim*EEGmerge.data)')')),'k','LineWidth',3);
p(4) = plot( max(abs(corr(EEG.srcpot',(Wfast*EEGmerge.data)')')),'b','LineWidth',3);
p(5) = plot( max(abs(corr(EEG.srcpot',(W1h*EEGmerge.data)')')),'r','LineWidth',3);
ylabel('Absolute Maximum Correlation','FontSize',16)
legend([p(5),p(4),p(3),p(2),p(1)],'PWC-ICA (1) Haar','FastICA','Ext. Infomax','Hilbert Analytic','Baseline')

%% Greedy
EEG.srcpot = data;
figure; hold on;
p(1) = errorbar([1:13],...
    mean(abs(corr(EEGmerge.srcpot',(EEGmerge.data)')')),...
    std(abs(corr(EEGmerge.srcpot',(EEGmerge.data)')')),...
    std(abs(corr(EEGmerge.srcpot',(EEGmerge.data)')')),'g','LineWidth',2);

p(2) = plot( greedyCorr(EEG.srcpot,WH*EEGmerge.data),'c','LineWidth',3);
p(3) = plot( greedyCorr(EEG.srcpot,Wim*EEGmerge.data),'k','LineWidth',3);
p(4) = plot( greedyCorr(EEG.srcpot,Wfast*EEGmerge.data),'b','LineWidth',3);
p(5) = plot( greedyCorr(EEG.srcpot,W1h*EEGmerge.data),'r','LineWidth',3);
ylabel('Absolute Maximum Correlation','FontSize',16)
xlabel('Source Number','FontSize',16)
legend([p(5),p(4),p(3),p(2),p(1)],'PWC-ICA (1) Haar','FastICA','Ext. Infomax','Hilbert Analytic','Baseline')

%% DIPOLE FITTING (Requires EEGLAB and DIPFIT extension)

key = {'Ext. Infomax',...
    'FastICA',...
    'PWC-ICA (1)',...
    'PWC-ICA (2)',...
    'PWC-ICA (4)',...
    'PWC-ICA (8)',...
    'PWC-ICA (1) Haar',...
    'PWC-ICA (2) Haar',...
    'PWC-ICA (4) Haar',...
    'PWC-ICA (8) Haar',...
    'PWC-ICA SG 3-7',...
    'PWC-ICA SG 2-5',...
    'PWC-ICA SG 5-9',...
    'PWC-ICA SG 7-21',...
    'PWC-ICA SG 11-51'...
    'Hilbert Analytic'
    }';

%% Infomax, FastICA
load(forwardParametersPath)
EEG.icaweights = Wim;
EEG.icawinv = inv(Wim);
EEG = pop_multifit(EEG, [1:64] ,'threshold',100,'plotopt',{'normlen' 'on'});

% Run once to warp position coordinates:
warp = traditionaldipfit(EEG.dipfit.coord_transform);
allposxyzsource = warp*[sourceCoords';ones(1,13)];
allposxyzsource = allposxyzsource(1:3,:);
allposxyzsource = allposxyzsource';
%

allRV(1,:) = cell2mat({ EEG.dipfit.model.rv });
for ii = 1:length(EEG.dipfit.model)
    allPosXYZ(:,ii,1) = EEG.dipfit.model(ii).posxyz(1,:)';
    allMomXYZ(:,ii,1) = EEG.dipfit.model(ii).momxyz(1,:)';
end
clear EEG;

load(forwardParametersPath)
EEG.icaweights = Wfast;
EEG.icawinv = inv(Wfast);
EEG = pop_multifit(EEG, [1:64] ,'threshold',100,'plotopt',{'normlen' 'on'});

allRV(2,:) = cell2mat({ EEG.dipfit.model.rv });
for ii = 1:length(EEG.dipfit.model)
    allPosXYZ(:,ii,2) = EEG.dipfit.model(ii).posxyz(1,:)';
    allMomXYZ(:,ii,2) = EEG.dipfit.model(ii).momxyz(1,:)';
end
clear EEG;

%% PWC-ICA
load(forwardParametersPath)
EEG.icaweights = W1;
EEG.icawinv = inv(W1);
EEG = pop_multifit(EEG, [1:64] ,'threshold',100,'plotopt',{'normlen' 'on'});


allRV(3,:) = cell2mat({ EEG.dipfit.model.rv });
for ii = 1:length(EEG.dipfit.model)
    allPosXYZ(:,ii,3) = EEG.dipfit.model(ii).posxyz(1,:)';
    allMomXYZ(:,ii,3) = EEG.dipfit.model(ii).momxyz(1,:)';
end
clear EEG;

load(forwardParametersPath)
EEG.icaweights = W2;
EEG.icawinv = inv(W2);
EEG = pop_multifit(EEG, [1:64] ,'threshold',100,'plotopt',{'normlen' 'on'});

allRV(4,:) = cell2mat({ EEG.dipfit.model.rv });
for ii = 1:length(EEG.dipfit.model)
    allPosXYZ(:,ii,4) = EEG.dipfit.model(ii).posxyz(1,:)';
    allMomXYZ(:,ii,4) = EEG.dipfit.model(ii).momxyz(1,:)';
end
clear EEG;

load(forwardParametersPath)
EEG.icaweights = W4;
EEG.icawinv = inv(W4);
EEG = pop_multifit(EEG, [1:64] ,'threshold',100,'plotopt',{'normlen' 'on'});

allRV(5,:) = cell2mat({ EEG.dipfit.model.rv });
for ii = 1:length(EEG.dipfit.model)
    allPosXYZ(:,ii,5) = EEG.dipfit.model(ii).posxyz(1,:)';
    allMomXYZ(:,ii,5) = EEG.dipfit.model(ii).momxyz(1,:)';
end
clear EEG;

load(forwardParametersPath)
EEG.icaweights = W8;
EEG.icawinv = inv(W8);
EEG = pop_multifit(EEG, [1:64] ,'threshold',100,'plotopt',{'normlen' 'on'});

allRV(6,:) = cell2mat({ EEG.dipfit.model.rv });
for ii = 1:length(EEG.dipfit.model)
    allPosXYZ(:,ii,6) = EEG.dipfit.model(ii).posxyz(1,:)';
    allMomXYZ(:,ii,6) = EEG.dipfit.model(ii).momxyz(1,:)';
end
clear EEG;

%% PWC-ICA Haar
load(forwardParametersPath)
EEG.icaweights = W1h;
EEG.icawinv = inv(W1h);
EEG = pop_multifit(EEG, [1:64] ,'threshold',100,'plotopt',{'normlen' 'on'});

allRV(7,:) = cell2mat({ EEG.dipfit.model.rv });
for ii = 1:length(EEG.dipfit.model)
    allPosXYZ(:,ii,7) = EEG.dipfit.model(ii).posxyz(1,:)';
    allMomXYZ(:,ii,7) = EEG.dipfit.model(ii).momxyz(1,:)';
end
clear EEG;

load(forwardParametersPath)
EEG.icaweights = W2h;
EEG.icawinv = inv(W2h);
EEG = pop_multifit(EEG, [1:64] ,'threshold',100,'plotopt',{'normlen' 'on'});

allRV(8,:) = cell2mat({ EEG.dipfit.model.rv });
for ii = 1:length(EEG.dipfit.model)
    allPosXYZ(:,ii,8) = EEG.dipfit.model(ii).posxyz(1,:)';
    allMomXYZ(:,ii,8) = EEG.dipfit.model(ii).momxyz(1,:)';
end
clear EEG;

load(forwardParametersPath)
EEG.icaweights = W4h;
EEG.icawinv = inv(W4h);
EEG = pop_multifit(EEG, [1:64] ,'threshold',100,'plotopt',{'normlen' 'on'});

allRV(9,:) = cell2mat({ EEG.dipfit.model.rv });
for ii = 1:length(EEG.dipfit.model)
    allPosXYZ(:,ii,9) = EEG.dipfit.model(ii).posxyz(1,:)';
    allMomXYZ(:,ii,9) = EEG.dipfit.model(ii).momxyz(1,:)';
end
clear EEG;

load(forwardParametersPath)
EEG.icaweights = W8h;
EEG.icawinv = inv(W8h);
EEG = pop_multifit(EEG, [1:64] ,'threshold',100,'plotopt',{'normlen' 'on'});

allRV(10,:) = cell2mat({ EEG.dipfit.model.rv });
for ii = 1:length(EEG.dipfit.model)
    allPosXYZ(:,ii,10) = EEG.dipfit.model(ii).posxyz(1,:)';
    allMomXYZ(:,ii,10) = EEG.dipfit.model(ii).momxyz(1,:)';
end
clear EEG;

%% Hilbert:
load(forwardParametersPath)
EEG.icaweights = WH;
EEG.icawinv = inv(WH);
EEG = pop_multifit(EEG, [1:64] ,'threshold',100,'plotopt',{'normlen' 'on'});

allRV(16,:) = cell2mat({ EEG.dipfit.model.rv });
for ii = 1:length(EEG.dipfit.model)
    allPosXYZ(:,ii,16) = EEG.dipfit.model(ii).posxyz(1,:)';
    allMomXYZ(:,ii,16) = EEG.dipfit.model(ii).momxyz(1,:)';
end
clear EEG;

%% Scatter plots for paper
key1 = key([1,2,16,7]);
for nn = 1:16
    clear dists corrs min1corr min2dist
    switch nn
        case 1; WW = Wim; case 2; WW = Wfast;
        case 3; WW = W1; case 4; WW = W2; case 5; WW = W4; case 6; WW = W8;
        case 7; WW = W1h; case 8; WW = W2h; case 9; WW = W4h; case 10; WW = W8h;
        case 11; WW = Wsg3_7; case 12; WW = Wsg2_5; case 13; WW = Wsg5_9; case 14; WW = Wsg7_21; case 15; WW = Wsg11_51;
        case 16; WW = WH;
    end
    for ii = 1:13
        dists(ii,:) = sqrt(sum(bsxfun(@minus,squeeze(allPosXYZ(:,:,nn))',squeeze(allposxyzsource(ii,:))).^2,2))';
    end
    corrs = abs(corr(data',(WW*EEGmerge.data)'));
    dists = dists(:,find(allRV(nn,:)<=.15));
    corrs = corrs(:,find(allRV(nn,:)<=.15));
    [~,min1ind] = min(dists);
    [~,min2ind] = max(corrs);
    for ii = 1:length(min1ind)
        min1corr(ii) = corrs(min1ind(ii),ii);
    end
    scatData1{nn} = [min1corr;min(dists)];
    for ii = 1:length(min2ind)
        min2dist(ii) = dists(min2ind(ii),ii);
    end
    scatData2{nn} = [max(corrs);min2dist];
    
end

figure; hold on; 
scatter(scatData1{1}(1,:),scatData1{1}(2,:),100,'kx','LineWidth',2.5);
scatter(scatData1{2}(1,:),scatData1{2}(2,:),100,'b+','LineWidth',2.5);
scatter(scatData1{16}(1,:),scatData1{16}(2,:),100,'cs','LineWidth',2.5);
scatter(scatData1{7}(1,:),scatData1{7}(2,:),100,'ro','LineWidth',2.5);
% scatter(scatData1{1}(1,:),scatData1{1}(2,:),(1-allRV(1,find(allRV(1,:)<=.15)))*100,'kx','LineWidth',2.5);
% scatter(scatData1{2}(1,:),scatData1{2}(2,:),(1-allRV(2,find(allRV(2,:)<=.15)))*100,'b+','LineWidth',2.5);
% scatter(scatData1{7}(1,:),scatData1{7}(2,:),(1-allRV(7,find(allRV(7,:)<=.15)))*100,'ro','LineWidth',2.5);
% scatter(scatData1{16}(1,:),scatData1{16}(2,:),(1-allRV(16,find(allRV(16,:)<=.15)))*100,'cs','LineWidth',2.5);
xlim([0,1]);
ylim([0,160]);
% title('Distance vs. Correlation of Component and Closest Source')
% xlabel('Correlation Factor between Component and its Closest Source')
% ylabel(sprintf('Euclidean Distance between\n Component and its Closest Source (mm)'))
xlabel('Correlation','FontSize',14)
ylabel('Euclidean Distance (mm)','FontSize',14)
set(gca,'FontSize',14)
legend(key1)

figure; hold on; 
scatter(scatData2{1}(1,:),scatData2{1}(2,:),100,'kx','LineWidth',2.5);
scatter(scatData2{2}(1,:),scatData2{2}(2,:),100,'b+','LineWidth',2.5);
scatter(scatData2{16}(1,:),scatData2{16}(2,:),100,'cs','LineWidth',2.5);
scatter(scatData2{7}(1,:),scatData2{7}(2,:),100,'ro','LineWidth',2.5);
% scatter(scatData2{1}(1,:),scatData2{1}(2,:),(1-allRV(1,find(allRV(1,:)<=.15)))*100,'kx','LineWidth',2.5);
% scatter(scatData2{2}(1,:),scatData2{2}(2,:),(1-allRV(2,find(allRV(2,:)<=.15)))*100,'b+','LineWidth',2.5);
% scatter(scatData2{7}(1,:),scatData2{7}(2,:),(1-allRV(7,find(allRV(7,:)<=.15)))*100,'ro','LineWidth',2.5);
% scatter(scatData2{16}(1,:),scatData2{16}(2,:),(1-allRV(16,find(allRV(16,:)<=.15)))*100,'cs','LineWidth',2.5);
xlim([0,1]);
ylim([0,160]);
% title(sprintf('Distance vs. Correlation of Components\n and Most Correlated Source'))
% xlabel(sprintf('Correlation Factor between Component\n and its Most Correlated Source'))
% ylabel(sprintf('Euclidean Distance between\n Component and its Most Correlated Source (mm)'))
xlabel('Correlation','FontSize',14)
ylabel('Euclidean Distance (mm)','FontSize',14)
set(gca,'FontSize',14)
legend(key1)

%% Do Comparisons (closest dipoles under 15% RV and most corellated dipoles under 15% RV)


for nn = [1:10,16]
    clear dists corrs min1corr min2dist
    switch nn
        case 1; WW = Wim; case 2; WW = Wfast;
        case 3; WW = W1; case 4; WW = W2; case 5; WW = W4; case 6; WW = W8;
        case 7; WW = W1h; case 8; WW = W2h; case 9; WW = W4h; case 10; WW = W8h;
        case 11; WW = Wsg3_7; case 12; WW = Wsg2_5; case 13; WW = Wsg5_9; case 14; WW = Wsg7_21; case 15; WW = Wsg11_51;
        case 16; WW = WH;
    end
    for ii = 1:13
        dists(ii,:) = sqrt(sum(bsxfun(@minus,squeeze(allPosXYZ(:,:,nn))',squeeze(allposxyzsource(ii,:))).^2,2))';
    end
    corrs = abs(corr(data',(WW*EEGmerge.data)'));
    dists = dists(:,find(allRV(nn,:)<=.15));
    corrs = corrs(:,find(allRV(nn,:)<=.15));
    [~,min1ind] = min(dists);
    [~,min2ind] = max(corrs);
    for ii = 1:length(min1ind)
        min1corr(ii) = corrs(min1ind(ii),ii);
    end
    scatData1{nn} = [min1corr;min(dists)];
    for ii = 1:length(min2ind)
        min2dist(ii) = dists(min2ind(ii),ii);
    end
    scatData2{nn} = [max(corrs);min2dist];
    
end

%%
key1 = key([1,2,16,7]);
for nn = 1:16
    clear dists corrs min1corr min2dist
    switch nn
        case 1; WW = Wim; case 2; WW = Wfast;
        case 3; WW = W1; case 4; WW = W2; case 5; WW = W4; case 6; WW = W8;
        case 7; WW = W1h; case 8; WW = W2h; case 9; WW = W4h; case 10; WW = W8h;
        case 11; WW = Wsg3_7; case 12; WW = Wsg2_5; case 13; WW = Wsg5_9; case 14; WW = Wsg7_21; case 15; WW = Wsg11_51;
        case 16; WW = WH;
    end
    for ii = 1:13
        dists(ii,:) = sqrt(sum(bsxfun(@minus,squeeze(allPosXYZ(:,:,nn))',squeeze(allposxyzsource(ii,:))).^2,2))';
    end
    corrs = abs(corr(data',(WW*EEGmerge.data)'));
    dists = dists(:,find(allRV(nn,:)<=1));
    corrs = corrs(:,find(allRV(nn,:)<=1));
    [~,min1ind] = min(dists);
    [~,min2ind] = max(corrs);
    for ii = 1:length(min1ind)
        min1corr(ii) = corrs(min1ind(ii),ii);
    end
    scatData1{nn} = [min1corr;min(dists)];
    for ii = 1:length(min2ind)
        min2dist(ii) = dists(min2ind(ii),ii);
    end
    scatData2{nn} = [max(corrs);min2dist];
    
end

figure; hold on; 
scatter(scatData1{1}(1,:),scatData1{1}(2,:),100,'kx','LineWidth',2.5);
scatter(scatData1{2}(1,:),scatData1{2}(2,:),100,'b+','LineWidth',2.5);
scatter(scatData1{16}(1,:),scatData1{16}(2,:),100,'cs','LineWidth',2.5);
scatter(scatData1{7}(1,:),scatData1{7}(2,:),100,'ro','LineWidth',2.5);
% scatter(scatData1{1}(1,:),scatData1{1}(2,:),(1-allRV(1,:))*100,'kx','LineWidth',2.5);
% scatter(scatData1{2}(1,:),scatData1{2}(2,:),(1-allRV(2,:))*100,'b+','LineWidth',2.5);
% scatter(scatData1{7}(1,:),scatData1{7}(2,:),(1-allRV(7,:))*100,'ro','LineWidth',2.5);
% scatter(scatData1{16}(1,:),scatData1{16}(2,:),(1-allRV(16,:))*100,'cs','LineWidth',2.5);
title('Distance vs. Correlation of Components and Closest Source')
% xlabel('Correlation Factor between Component and its Closest Source')
% ylabel(sprintf('Euclidean Distance between\n Component and its Closest Source (mm)'))
xlabel('Correlation','FontSize',14)
ylabel('Euclidean Distance (mm)','FontSize',14)
set(gca,'FontSize',14)
xlim([0,1]);
ylim([0,160])
legend(key1)

figure; hold on; 
scatter(scatData2{1}(1,:),scatData2{1}(2,:),100,'kx','LineWidth',2.5);
scatter(scatData2{2}(1,:),scatData2{2}(2,:),100,'b+','LineWidth',2.5);
scatter(scatData2{16}(1,:),scatData2{16}(2,:),100,'cs','LineWidth',2.5);
scatter(scatData2{7}(1,:),scatData2{7}(2,:),100,'ro','LineWidth',2.5);
% scatter(scatData2{1}(1,:),scatData2{1}(2,:),(1-allRV(1,:))*100,'kx');
% scatter(scatData2{2}(1,:),scatData2{2}(2,:),(1-allRV(2,:))*100,'b+');
% scatter(scatData2{7}(1,:),scatData2{7}(2,:),(1-allRV(7,:))*100,'ro');
% scatter(scatData2{16}(1,:),scatData2{16}(2,:),(1-allRV(16,:))*100,'cs');
title(sprintf('Distance vs. Correlation of Components\n and Most Correlated Source'))
% xlabel(sprintf('Correlation Factor between Component\n and its Most Correlated Source'))
% ylabel(sprintf('Euclidean Distance between\n Component and its Most Correlated Source (mm)'))
xlabel('Correlation','FontSize',14)
ylabel('Euclidean Distance (mm)','FontSize',14)
set(gca,'FontSize',14)
xlim([0,1]);
ylim([0,160])
legend(key1)