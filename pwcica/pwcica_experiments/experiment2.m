%PWC-ICA Paper; 
% Experiment 2;
% Instructions:
% Either load pwcica_exp2.mat, or if you have the SIFT toolbox you may run
% the code in the first section to generate a randomized version of the
% data.

load('pwcica_exp2.mat');
%% (OPTIONAL, need SIFT) Generate simulated sources.
% Example 2:  System of coupled oscillators 
% from (Ex 3.1, eq 11-15) Schelter B, Timmer J, Eichler M (2009) Assessing the strength of directed 
% influences among neural signals using renormalized partial directed coherence. 
% Journal of neuroscience methods 179:121-30 Available at: http://www.ncbi.nlm.nih.gov/pubmed/19428518.

SamplingRate = 200;      % Hz

Nl = 500;                % length of each epoch (samples)
Nr = 100;                % number of trials (realisations)
ndisc = 1000;            % number of samples to discard from VAR model (startup transients)
ModelOrder = 3;          % VAR model order
sourceCount = 5;

% write the system of equations for coupled damped oscillators
expr = { ...
    'x1(t) = 0.9*x1(t-1)  + 0.3*x2(t-2)  + e1(t)' ...
    'x2(t) = 1.3*x2(t-1)  + -0.8*x2(t-2) + e2(t)' ...
    'x3(t) = 0.3*x1(t-2)  + 0.6*x2(t-1)  + e3(t)' ...
    'x4(t) = -0.7*x4(t-3) + -0.7*x1(t-3) + 0.3*x5(t-3) + e4(t)' ...
    'x5(t) = 1*x5(t-1)    + -0.4*x5(t-2) + 0.3*x4(t-2) + e5(t)' ...
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

%% Setup Mixing Experiment (Optional: pwcica_exp1.mat includes the randomized mixing matrices used for paper. Run to generate new experiment).

% % Zero mean the source data.
% data = bsxfun(@minus,data,mean(data,2));

% Define random-mixings
% if ~exist('mixes')
%     mixes = cell(10,1);
    for ii = 1:20
        mixes{ii} = randn(sourceCount);%/sqrt(10);
        while cond(mixes{ii}) > 2*sourceCount;
        mixes{ii} = randn(sourceCount);%*sqrt(10);
        end
    end

% end

% Build signal datasets:
sig = cell(size(mixes));
for ii = 1:length(sig)
    sig{ii} = reshape(mixes{ii}*reshape(data,size(data,1),size(data,2)*size(data,3)),size(data));
    sig{ii} = bsxfun(@minus,sig{ii},mean(sig{ii},2));
end

%% Run extended Infomax ICA (requires runica() found in EEGLAB)
if ~exist('Wim')
    Wim = cell(size(mixes));
end

for ii = 1:20
    
    tic;
    [Wim{ii},Sim] = runica(reshape(sig{ii},size(data,1),size(data,2)*size(data,3)),'extended',1);
    timeElapsed.im(ii) = toc;
    Wim{ii} = Wim{ii}*Sim;
    
end

%% Run fastica (requires fastica(), found in fieldtrip plugin for EEGLAB)

if ~exist('Wfast')
    Wfast = cell(size(mixes));
end

for ii = 1:20
    
    tic;
    [~,Wfast{ii}] = fastica(reshape(sig{ii},size(data,1),size(data,2)*size(data,3)),'approach','symm');
    timeElapsed.fast(ii) = toc;
    
end

%% Run pwcica
if ~exist('W1')
    W1 = cell(size(mixes));
    W2 = cell(size(mixes));
    W4 = cell(size(mixes));
    W8 = cell(size(mixes));
    metas1 = cell(size(mixes));
    metas2 = cell(size(mixes));
    metas4 = cell(size(mixes));
    metas8 = cell(size(mixes));
end

for ii = 1:20
    tic;
    [W1{ii},metas1{ii}] = pwcica(sig{ii},'TimeInvariant',0,'SamplingRate',SamplingRate,'Level',1,'NoPaths',1);
    timeElapsed.pwcica1(ii) = toc;
    tic;
    [W2{ii},metas2{ii}] = pwcica(sig{ii},'TimeInvariant',0,'SamplingRate',SamplingRate,'Level',2,'NoPaths',1);
    timeElapsed.pwcica2(ii) = toc;
    tic;
    [W4{ii},metas4{ii}] = pwcica(sig{ii},'TimeInvariant',0,'SamplingRate',SamplingRate,'Level',4,'NoPaths',1);
    timeElapsed.pwcica4(ii) = toc;
    tic;
    [W8{ii},metas8{ii}] = pwcica(sig{ii},'TimeInvariant',0,'SamplingRate',SamplingRate,'Level',8,'NoPaths',1);
    timeElapsed.pwcica8(ii) = toc;
end

%% Run pwcica using Haar wavelets.
if ~exist('W1h')
    W1h = cell(size(mixes));
    W2h = cell(size(mixes));
    W4h = cell(size(mixes));
    W8h = cell(size(mixes));
    metas1h = cell(size(mixes));
    metas2h = cell(size(mixes));
    metas4h = cell(size(mixes));
    metas8h = cell(size(mixes));
end

for ii = 1:20
    tic;
    [W1h{ii},metas1h{ii}] = pwcica(sig{ii},'TimeInvariant',1,'SamplingRate',SamplingRate,'Level',1,'NoPaths',1);
    timeElapsed.pwcica1h(ii) = toc;
    tic;
    [W2h{ii},metas2h{ii}] = pwcica(sig{ii},'TimeInvariant',1,'SamplingRate',SamplingRate,'Level',2,'NoPaths',1);
    timeElapsed.pwcica2h(ii) = toc;
    tic;
    [W4h{ii},metas4h{ii}] = pwcica(sig{ii},'TimeInvariant',1,'SamplingRate',SamplingRate,'Level',4,'NoPaths',1);
    timeElapsed.pwcica4h(ii) = toc;
    tic;
    [W8h{ii},metas8h{ii}] = pwcica(sig{ii},'TimeInvariant',1,'SamplingRate',SamplingRate,'Level',8,'NoPaths',1);
    timeElapsed.pwcica8h(ii) = toc;
end

%% Run pwcica using Hilbert transform...
if ~exist('WH')
    WH = cell(size(mixes));
    metasH = cell(size(mixes));
end

for ii = 1:20
    tic;
    [WH{ii},metasH{ii}] = pwcica(sig{ii},'TimeInvariant',1,'SamplingRate',SamplingRate,'Level',0,'NoPaths',1);
    timeElapsed.pwcicaH(ii) = toc;
end

%% New Amari Index calcs:
aI = zeros(15,20);
tE = zeros(15,20);
aIs= zeros(15,20);

% infomax
index = 1;
for jj = 1:20
    scaling = repmat(sqrt(mean(inv(Wim{jj}).^2))',[1 size(Wim{jj},1)]);
    aIs(index,jj) = amariIndex(scaling.*Wim{jj},inv(mixes{jj}));
    aI(index,jj) = amariIndex(Wim{jj},inv(mixes{jj}));
    tE(index,jj) = timeElapsed.im(jj);
end
% fastICA
index = 2;
for jj = 1:20
    scaling = repmat(sqrt(mean(inv(Wfast{jj}).^2))',[1 size(Wfast{jj},1)]);
    aIs(index,jj) = amariIndex(scaling.*Wfast{jj},inv(mixes{jj}));
    aI(index,jj) = amariIndex(Wfast{jj},inv(mixes{jj}));
    tE(index,jj) = timeElapsed.fast(jj);
end
% pwc-ica
for ii = 1:4
    index = 2+ii;
    varname = sprintf('%i',2^(ii-1));
    W = eval(['W',varname]);
    times = eval(['timeElapsed.pwcica',varname]);
    for jj = 1:20
        scaling = repmat(sqrt(mean(inv(W{jj}).^2))',[1 size(W{jj},1)]);
        aIs(index,jj) = amariIndex(scaling.*W{jj},inv(mixes{jj}));
        aI(index,jj) = amariIndex(W{jj},inv(mixes{jj}));
        tE(index,jj) = times(jj);
    end
    clear W times
end
% pwc-ica Haar
for ii = 1:4
    index = 6+ii;
    varname = sprintf('%ih',2^(ii-1));
    W = eval(['W',varname]);
    times = eval(['timeElapsed.pwcica',varname]);
    for jj = 1:20
        scaling = repmat(sqrt(mean(inv(W{jj}).^2))',[1 size(W{jj},1)]);
        aIs(index,jj) = amariIndex(scaling.*W{jj},inv(mixes{jj}));
        aI(index,jj) = amariIndex(W{jj},inv(mixes{jj}));
        tE(index,jj) = times(jj);
    end
    clear W times
end
% fastICA with Hilbert transform.
index = 11;
for jj = 1:20
    scaling = repmat(sqrt(mean(inv(WH{jj}).^2))',[1 size(WH{jj},1)]);
    aIs(index,jj) = amariIndex(scaling.*WH{jj},inv(mixes{jj}));
    aI(index,jj) = amariIndex(WH{jj},inv(mixes{jj}));
    tE(index,jj) = timeElapsed.pwcicaH(jj);
end

%% means and stds of Amari Indices:
mean(aIs')';
std(aIs')';