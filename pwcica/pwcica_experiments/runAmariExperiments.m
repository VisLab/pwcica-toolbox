%PWC-ICA Paper - this script generates the data for Table 1.
% It assumes you have already generated the random mixtures or are
% loading the mixtures used in the paper. The script produces the
% values and saves the values in a report file.
% 

%% Run experiment 1
modelFile = 'D:\Research\Amica\pwcica-toolbox\pwcica\pwcica_experiments\data\pwcica_exp1.mat';
reportFile = 'D:\Research\Amica\pwcica_exp1_report.mat'; 
if ~isempty(modelFile)
   load(modelFile);
   model = modelExp1;
else
   model = generateAmariExperiment(1);
end
experiment = doAmariExperiment(model, 'Uncoupled oscillator model');  %#ok<NASGU>
save(reportFile, 'experiment', '-v7.3');

%% Run experiment 2
modelFile = 'D:\Research\Amica\pwcica-toolbox\pwcica\pwcica_experiments\data\pwcica_exp2.mat';
reportFile = 'D:\Research\Amica\pwcica_exp2_report.mat'; 
if ~isempty(modelFile)
   load(modelFile);
   model = modelExp2;
else
   model = generateAmariExperiment(2);
end
experiment = doAmariExperiment(model, 'Uncoupled oscillator model');  %#ok<NASGU>
save(reportFile, 'experiment', '-v7.3');

%% Run experiment 3
modelFile = 'D:\Research\Amica\pwcica-toolbox\pwcica\pwcica_experiments\data\pwcica_exp3.mat';
reportFile = 'D:\Research\Amica\pwcica_exp3_report.mat'; 
if ~isempty(modelFile)
   load(modelFile);
   model = modelExp3;
else
   model = generateAmariExperiment(3);
end
experiment = doAmariExperiment(model, 'Dynamic oscillator model'); 
save(reportFile, 'experiment', '-v7.3');