%% Find the closest and most correlated sources (Experiment 5)
% 
% If you have not calculated pwcica_exp5.mat by running the correlation
% experiment, you must do so before starting this.

%% Initialize
dataDir = 'D:\Research\Amica\pwcica-toolbox\pwcica\pwcica_experiments\data'; 
forwardParamFile = 'D:\Research\Amica\pwcica-toolbox\pwcica\pwcica_experiments\data\forwardParameters';
inputFile = 'pwcica_exp5.mat';
experimentName = 'exp5';

%% Load the correlation test results
load([dataDir filesep inputFile]);

%% DIPOLE FITTING (Requires EEGLAB and DIPFIT extension)
allRV = cell(length(ICAMethods), 1);
allPosXYZ = cell(length(ICAMethods), 1);
allMomXYZ = cell(length(ICAMethods), 1);
warp = [];
for k = 1:length(ICAMethods);
    clear EEG;
    load(forwardParamFile)   % Load the source coordinates
    EEG.icaweights = W{k};
    EEG.icawinv = inv(W{k});
    EEG = pop_multifit(EEG, 1:64 , 'threshold', 100, ...
        'plotopt',{'normlen' 'on'}); 
    if isempty(warp)  % Set the warp if empty
        warp = traditionaldipfit(EEG.dipfit.coord_transform);
    end
    allRV{k} = cell2mat({EEG.dipfit.model.rv});
    model = EEG.dipfit.model;
    for ii = 1:length(model)
        allPosXYZ{k}(:, ii) = model(ii).posxyz(1, :)';
        allMomXYZ{k}(:, ii) = model(ii).momxyz(1, :)';
    end
end

%% Find the warp position coordinates for the actual sources
load(forwardParamFile)
allposxyzsource = warp*[sourceCoords'; ones(1, 13)];
allposxyzsource = allposxyzsource(1:3, :);
allposxyzsource = allposxyzsource';

%% Scatter plots for paper
eps = 0.15;
scatterTypes = {'Smallest distance', 'Smallest correlation', ...
         'Smallest distance, RV < eps', 'Smallest correlation, RV < eps'};
scatterData = cell(length(ICAMethods), length(scatterTypes)); 

data = EEGmerge.etc.srcPotentials;
for k = 1:length(ICAMethods);
    allPos = allPosXYZ{k};
    theDist = pdist2(allPos', allposxyzsource);
    theCorrs = abs(corr(data', (W{k}*EEGmerge.data)'))';

    [minDist, minind] = min(theDist, [], 2);
    [maxCorr, maxind] = max(theCorrs, [], 2);
    smallDist = theDist(allRV{k} <= eps, :);
    smallCorrs = theCorrs(allRV{k} <= eps, :);
    [minSmallDist, minSmallInd] = min(smallDist, [], 2);
    [maxSmallCorrs, maxSmallInd] = max(smallCorrs, [], 2);
    
    d1 = [minDist, minDist];
    d2 = [maxCorr, maxCorr];
    for j = 1:size(theDist, 1)
        d1(j, 1) = theCorrs(j, minind(j));
        d2(j, 2) = theDist(j, maxind(j));
    end
    d3 = [minSmallDist, minSmallDist];
    d4 = [maxSmallCorrs, maxSmallCorrs];
    for j = 1:size(smallDist, 1)
        d3(j, 1) = smallCorrs(j, minSmallInd(j));
        d4(j, 2) = smallDist(j, maxSmallInd(j));
    end
    scatterData{k, 1} = d1;
    scatterData{k, 2} = d2; 
    scatterData{k, 3} = d3;
    scatterData{k, 4} = d4;     
end

%% Set up the parameters for the plots  (add your method, another color and shape)
methods = {'Hilbert', 'Fastica', 'Infomax', 'Amica', 'Pwcica1h'};
mLegends = {'PWC-ICA(1) Haar', 'Amica', 'Ext Infomax', 'Fastica', 'Hilbert'}; 
mColors = [0, 1, 1; ...      % Hilbert
           0, 0.5, 0.5; ...  % Fastica
           0, 0, 0; ...      % Infomax
           0, 0, 1; ...      % Amica
           1, 0, 0];         % Pwcica1h
mMarkers = {'s' ; ...   % Hilbert
            'd'; ...  % Fastica
            'x'; ...  % Infomax
            '+'; ...  % Amica
            'o'};     % Pwcica1h
        
%% Create the plots
for stype = 1:size(scatterData, 2);
    figure('Color', [1, 1, 1], 'Name', scatterTypes{stype});  
    hold on;
    p = [];
    set(gca,'FontSize', 14);
    for k = 1:length(methods);
       m = find(strcmpi(ICAMethods, methods(k)));
       theData = scatterData{m, stype};
       q = plot(theData(:, 1), theData(:, 2), ...
                   'Marker', mMarkers{k}, 'MarkerEdgeColor', mColors(k, :), ...
                   'MarkerFaceColor', 'none', 'MarkerSize', 10, ...
                   'LineWidth', 2.5, 'LineStyle', 'none');
       p = [q p]; %#ok<AGROW>
    end
    xlim([0,1]);
    ylim([0,160]);
    xlabel('Correlation', 'FontSize', 14)
    ylabel('Euclidean Distance (mm)', 'FontSize', 14)
    legend(p, mLegends)
    hold off
end
