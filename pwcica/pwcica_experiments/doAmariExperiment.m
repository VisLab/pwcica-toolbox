function experiment = doAmariExperiment(model, experimentName)
%PWC-ICA Paper - this script generates the data for Table 1.
% It assumes you have already generated the random mixtures or are
% loading the mixtures used in the paper. The script produces the
% values and saves the values in a report file.
% 
% You can generate the experiments by running:
%  runAmari

%% Initialize the script - substitute your own mixture files if desired.
experiment = struct('name', experimentName, 'mixes',[], ...
                    'ICAMethods', [], 'Ws', [], ...
                    'amari', [], 'aMean', [], 'aStd', [], ...
                    'amariScaled', [], 'asMean', [], 'asStd', [], ...
                    'pValues', [], 'methodTimes', []);
ICAMethods = {};  

%% Initialize the model information
mixes = model.mixes;
numTrials = length(mixes);
data = model.data;
samplingRate = model.samplingRate;
%% Build signal datasets 
sig = cell(size(mixes));
for ii = 1:length(sig)
    sig{ii} = reshape(mixes{ii}*reshape(data, size(data,1), ...
                     size(data,2)*size(data,3)), size(data));
    sig{ii} = bsxfun(@minus,sig{ii}, mean(sig{ii}, 2));
end

%% Run extended Infomax ICA (requires runica() found in EEGLAB)
if ~exist('WInfomax', 'var')
    WInfomax = cell(size(mixes));
    ICAMethods{end+1} = 'Infomax';
end

for ii = 1:numTrials
    tic;
    [WInfomax{ii}, SInfomax] = runica(reshape(sig{ii},size(data,1),size(data,2)*size(data,3)),'extended',1);
    timeElapsed.Infomax(ii) = toc;
    WInfomax{ii} = WInfomax{ii}*SInfomax;  
end

%% Run amica (requires runamica15 found in the EEGLAB amica plugin)
if ~exist('WAmica', 'var')
    WAmica = cell(size(mixes));
    ICAMethods{end+1} = 'Amica';
end

for ii = 1:numTrials
    tic;
    [WAmica{ii}, SAmica] = runamica15( ...
        reshape(sig{ii}, size(data, 1), size(data, 2)*size(data, 3)), ...
        'outdir', 'amicaExp1');
    timeElapsed.Amica(ii) = toc;
    WAmica{ii} = WAmica{ii}*SAmica;  
end

%% Run fastica (requires fastica(), found in fieldtrip plugin for EEGLAB)
if ~exist('WFastica', 'var')
    WFastica = cell(size(mixes));
    ICAMethods{end+1} = 'Fastica';
end

for ii = 1:numTrials  
    tic;
    [~, WFastica{ii}] = fastica(reshape(sig{ii}, size(data, 1), ...
                   size(data, 2)*size(data, 3)),'approach', 'symm');
    timeElapsed.Fastica(ii) = toc;   
end

%% Run pwcica
if ~exist('WPwcica1', 'var')
    WPwcica1 = cell(size(mixes));
    WPwcica2 = cell(size(mixes));
    WPwcica4 = cell(size(mixes));
    WPwcica8 = cell(size(mixes));
    for k = [1, 2, 4, 8]
      ICAMethods{end + 1} = ['Pwcica' num2str(k)];   %#ok<AGROW>
    end
end

for ii = 1:numTrials
    tic;
    WPwcica1{ii} = pwcica(sig{ii}, 'TimeInvariant', 0, ...
                  'samplingRate', samplingRate, 'Level', 1, 'NoPaths', 1);
    timeElapsed.Pwcica1(ii) = toc;
    tic;
    WPwcica2{ii} = pwcica(sig{ii}, 'TimeInvariant', 0, ...
                    'samplingRate', samplingRate, 'Level', 2, 'NoPaths', 1);
    timeElapsed.Pwcica2(ii) = toc;
    tic;
    WPwcica4{ii} = pwcica(sig{ii}, 'TimeInvariant', 0, ...
                    'samplingRate', samplingRate, 'Level', 4, 'NoPaths', 1);
    timeElapsed.Pwcica4(ii) = toc;
    tic;
    WPwcica8{ii} = pwcica(sig{ii}, 'TimeInvariant', 0, ...
                   'samplingRate', samplingRate, 'Level', 8, 'NoPaths', 1);
    timeElapsed.Pwcica8(ii) = toc;
end

%% Run pwcica using Haar wavelets.
if ~exist('WPwcica1h', 'var')
    WPwcica1h = cell(size(mixes));
    WPwcica2h = cell(size(mixes));
    WPwcica4h = cell(size(mixes));
    WPwcica8h = cell(size(mixes));
    for k = [1, 2, 4, 8]
      ICAMethods{end + 1} = ['Pwcica' num2str(k) 'h'];   %#ok<AGROW>
    end
end

for ii = 1:numTrials
    tic;
    WPwcica1h{ii} = pwcica(sig{ii}, 'TimeInvariant', 1, ...
                      'samplingRate', samplingRate, 'Level', 1, 'NoPaths', 1);
    timeElapsed.Pwcica1h(ii) = toc;
    tic;
    WPwcica2h{ii} = pwcica(sig{ii}, 'TimeInvariant', 1, ...
                      'samplingRate', samplingRate, 'Level', 2, 'NoPaths', 1);
    timeElapsed.Pwcica2h(ii) = toc;
    tic;
    WPwcica4h{ii} = pwcica(sig{ii},'TimeInvariant', 1, ...
                     'samplingRate', samplingRate, 'Level', 4, 'NoPaths', 1);
    timeElapsed.Pwcica4h(ii) = toc;
    tic;
    WPwcica8h{ii} = pwcica(sig{ii}, 'TimeInvariant', 1, ...
                     'samplingRate', samplingRate, 'Level', 8, 'NoPaths', 1);
    timeElapsed.Pwcica8h(ii) = toc;
end

%% Run pwcica using Hilbert transform...
if ~exist('WHilbert', 'var')
    WHilbert = cell(size(mixes));
    ICAMethods{end+1} = 'Hilbert';
end

for ii = 1:numTrials
    tic;
    WHilbert{ii} = pwcica(sig{ii}, 'TimeInvariant', 1, ...
                   'samplingRate', samplingRate, 'Level', 0, 'NoPaths', 1);
    timeElapsed.Hilbert(ii) = toc;
end

%% New Amari Index calculations
numMethods = length(ICAMethods);
amari = zeros(numMethods, numTrials);
methodTimes = zeros(numMethods, numTrials);
amariScaled = zeros(numMethods, numTrials);

Ws = cell(numMethods, 1);
for k = 1:numMethods
    W = eval(['W', ICAMethods{k}]);
    Ws{k} = W;
    for jj = 1:numTrials
        scaling = repmat(sqrt(mean(inv(W{jj}).^2))', [1 size(W{jj}, 1)]);
        amariScaled(k, jj) = amariIndex(scaling.*W{jj}, inv(mixes{jj}));
        amari(k, jj) = amariIndex(W{jj}, inv(mixes{jj}));
        methodTimes(k, jj) = timeElapsed.(ICAMethods{k})(jj);
    end
end

%% Output the table of times
aMean = mean(amari, 2);
aStd = std(amari, 0, 2);
asMean = mean(amariScaled, 2);
asStd = std(amariScaled, 0, 2);
fprintf('Method    Mean     Std   ScaledMean  ScaledStD\n');
for k = 1:numMethods
    fprintf('%s\t  %5.3f\t  %5.3f  %5.3f\t  %5.3f\n',  ICAMethods{k}, ...
        aMean(k), aStd(k), asMean(k), asStd(k));
end

%% Calculate the table of p-Values using one-sided paired t-test
significance = 0.05;
pValues = ones(numMethods, numMethods);
for k = 1:numMethods
    for j = 1:numMethods
       [~, pValues(k, j)] = ttest(amari(k, :)', amari(j, :)', ...
                                  'Tail', 'left', 'Alpha', significance);     
    end
end

%% Output the table of p-Values 
for k = 1:numMethods
    fprintf('\nMethods better than %s\n', ICAMethods{k});
    for j = 1:numMethods
        if ~isnan(pValues(j, k)) && pValues(j, k) < significance
            fprintf('   %s\t %g\n',  ICAMethods{j}, pValues(j, k));
        end
    end
end

%% Return the output in a structure
experiment.mixes = mixes;
experiment.ICAMethods = ICAMethods;
experiment.Ws = Ws;
experiment.amari = amari;
experiment.aMean = aMean;
experiment.aStd = aStd;
experiment.amariScaled = amariScaled;
experiment.asMean = asMean;
experiment.asStd = asStd;
experiment.pValues = pValues;
experiment.methodTimes = methodTimes;