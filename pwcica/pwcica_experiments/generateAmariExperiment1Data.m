%% Generate simulated sources for 10-variate un-coupled damped oscillators
% Example 1a: 10-variate un-coupled damped oscillators
% Adapted from SIFT, by Tim Mullen. This program needs sift to run


SamplingRate = 200;      % Hz
Nl = 500;                % length of each epoch (samples)
Nr = 100;                % number of trials (realisations)
ndisc = 1000;            % number of samples to discard from VAR model (startup transients)
ModelOrder = 2;          % VAR model order
% ModelOrder = 4;        % Model Order we will use
f0 = 10;                 % central oscillation frequency
fr = [2,3,5,7,11,13,17,19,23,29];
tau = 5+randn(1,10)*3;
ind = find(tau < 0); 
tau(ind) = 2+3*abs(randn(1));
sourceCount = 10;
numberMixes = 20;
expr = {...
    ['x1(t) = ' sim_dampedOscillator(fr(1),tau(1),SamplingRate,1) ' + e1(t)'] ... 
    ['x2(t) = ' sim_dampedOscillator(fr(2),tau(2),SamplingRate,2) ' + e2(t)'] ...
    ['x3(t) = ' sim_dampedOscillator(fr(3),tau(3),SamplingRate,3) ' + e3(t)'] ...
    ['x4(t) = ' sim_dampedOscillator(fr(4),tau(4),SamplingRate,4) ' + e4(t)'] ...
    ['x5(t) = ' sim_dampedOscillator(fr(5),tau(5),SamplingRate,5) ' + e5(t)'] ...
    ['x6(t) = ' sim_dampedOscillator(fr(6),tau(6),SamplingRate,6) ' + e6(t)'] ...
    ['x7(t) = ' sim_dampedOscillator(fr(7),tau(7),SamplingRate,7) ' + e7(t)'] ...
    ['x8(t) = ' sim_dampedOscillator(fr(8),tau(8),SamplingRate,8) ' + e8(t)'] ...
    ['x9(t) = ' sim_dampedOscillator(fr(9),tau(9),SamplingRate,9) ' + e9(t)'] ...
    ['x10(t) = ' sim_dampedOscillator(fr(10),tau(10),SamplingRate,10) ' + e10(t)'] ...
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
mixes = cell(numberMixes, 1);
for ii = 1:numberMixes
    mixes{ii} = randn(sourceCount); %/sqrt(10);
    while cond(mixes{ii}) > 2*sourceCount;
        mixes{ii} = randn(sourceCount);%*sqrt(10);
    end
end

% end