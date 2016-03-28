%% This script runs the experiment for Fig. 4 of the paper 

%% Initialize Initialize the script
EEGFile = 'pwcica_exp4.mat';
outputDir = 'D:\Research\Amica\data'; 
resultsFile = 'CorrelationTest.mat';
experimentName = 'exp4a';

%% Initialize the values
ICAMethods = {};
W = {};

%% Load the data
if ~isempty(EEGFile)
   load(EEGFile);
else
   EEGmerge = generateForwardExperiment();
end

%% Run extended Infomax ICA (requires runica() found in EEGLAB)
tic;
ICAMethods{end+1} = 'Infomax';
EEG = pop_runica(EEGmerge, 'extended', 1);
W{end + 1} = EEG.icaweights*EEG.icasphere;
timeElapsed.Infomax = toc;

%% Run amica (requires runamica15 found in the EEGLAB amica plugin)
tic;
ICAMethods{end+1} = 'Amica';
[WAmica, SAmica] = runamica15(EEGmerge.data, ...
           'outdir', [outputDir filesep 'amica' experimentName]);
W{end + 1} = WAmica*SAmica;
timeElapsed.Amica = toc;

%% Run FastICA (requires fastica(), found in fieldtrip plugin for EEGLAB)
tic;
ICAMethods{end+1} = 'Fastica';
[~, W{end + 1}] = fastica(EEGmerge.data, 'approach', 'symm');
timeElapsed.Fastica = toc;

%% Run PWCICA
WPwcicaOrders = [1, 2, 4, 8, 64];
for k = WPwcicaOrders
   tic;
   methodName = ['Pwcica' num2str(k)];
   ICAMethods{end + 1} = methodName;  %#ok<SAGROW>
   W{end + 1} = pwcica(EEGmerge, 'Level', k, ...
                      'TimeInvariant', 0, 'NoPaths', 1); %#ok<SAGROW>
   timeElapsed.(methodName) = toc;
end

%% Run PWCICA with Haar
WPwcicahOrders = [1, 2, 4, 8, 64];
for k = WPwcicahOrders
   tic;
   methodName = ['Pwcica' num2str(k) 'h'];
   ICAMethods{end + 1} = methodName;  %#ok<SAGROW>
   W{end + 1} = pwcica(EEGmerge, 'Level', k, ...
                      'TimeInvariant', 1, 'NoPaths', 1); %#ok<SAGROW>
   timeElapsed.(methodName) = toc;
end

%% Run Hilbert transform approach
tic;
ICAMethods{end + 1} = 'Hilbert';
W{end + 1} = pwcica(EEGmerge, 'Level', 0, ...
                    'TimeInvariant', 1, 'NoPaths', 1);
timeElapsed.Hilbert = toc;

%% Specification of what methods to plot in what colors
methods = {'Hilbert', 'Fastica', 'Infomax', 'Amica', 'Pwcica1h'};
mLegends = {'PWC-ICA(1) Haar', 'Amica', 'Ext Infomax', 'Fastica', 'Hilbert', 'Baseline'}; 
mColors = [0, 1, 1; ...     % Hilbert
           0, 0.5, 0.5; ... % Fastica
           0, 0, 0; ...     % Infomax
           0, 0, 1; ...     % Amica
           1, 0, 0];        % Pwcica1h

%% Non-greedy plot     
figure 
hold on;
clear p;
p(1) = errorbar(1:13, mean(abs(corr(data',(EEGmerge.data)')')),...
               std(abs(corr(data',(EEGmerge.data)')')), ...
               'Color', [0.8, 0.8, 0.8], 'LineWidth', 2);
for k = 1:length(methods);
   m = find(strcmpi(ICAMethods, methods(k)));
   q = plot( max(abs(corr(data',(W{m}*EEGmerge.data)')')), ...
                     'Color', mColors(k, :), 'LineWidth',3); 
   p = [q p]; %#ok<AGROW>
end
ylabel('Absolute Maximum Correlation','FontSize',16)
legend(p, mLegends)
hold off

%% Greedy plot (add your method and a color)
figure 
hold on;
clear p;
p(1) = errorbar(1:13, mean(abs(corr(data',(EEGmerge.data)')')),...
               std(abs(corr(data',(EEGmerge.data)')')), ...
               'Color', [0.8, 0.8, 0.8], 'LineWidth', 2);
for k = 1:length(methods);
   m = find(strcmpi(ICAMethods, methods(k)));
   q = plot(greedyCorr(data, W{m}*EEGmerge.data), ...
                     'Color', mColors(k, :), 'LineWidth',3); 
   p = [q p]; %#ok<AGROW>

end
ylabel('Absolute Maximum Correlation','FontSize',16)
xlabel('Source Number','FontSize',16)
legend(p, mLegends)
hold off

%% Save the results
correlationResults.data = data;
correlationResults.EEGmerge = EEGmerge;
correlationResults.ICAMethods = ICAMethods;
correlationResults.W = W;
correlationResults.timeElapsed = timeElapsed;
save([outputDir filesep resultsFile], 'correlationResults', '-v7.3');
 
 