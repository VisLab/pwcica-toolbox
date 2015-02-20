% pwcica -- Pair-wise Complex ICA. 
%   Transforms time series data into a complex phase space, performs one of
%   several complex ICA transforms on SEQUENTIALLY ORDERED DATA, then 
%   projects complex demixing matrix back to base space resulting in a 
%   real-valued demixing matrix.
%
%   INPUTS: 
%     data - 
%     1) A [n,N,K] real-valued array of epoched time series data. n is
%        the dimension of individual observations. N is the number of
%        observations in an epoch. K is the total number of epochs.
% 
%     2) A [n,N] real-valued array of time series data. n is the dimension
%        of individual observations. N is the total number of observations.
% 
%     3) An EEGLAB EEG structure. In this case, pwcica2 assumes that data
%        is an [n,N,K] array contained in EEG.data.
%     
%     varargin - 
%     Optional key value pairs:
%     ** 'key' , ( [default] , option , option )
%     
%     -> 'ComplexICAMethod' , ( ['fastICA'] 'ebm' robust' @customfunction )
%           The complex ICA method to be used in complex phase space. NOTE:
%           use 'infomax' only with 'Real' set to 1.
%           ***!!! SEE NOTE BELOW ON COMPLEX ICA PACKAGES !!!***
%           ALTERNATELY: the value of ComplexICAMethod may be specified as a
%           function handle for a custom complex ICA method satisfying:
%           >>    complexW = @customfunction(complexData);
%           where complexData is the [n,N] array of complex signal data,
%           and complexW is the output [n,n] complex demixing matrix. Note,
%           pwcica automatically (pseudo) whitens complexData input.
%           EXAMPLE: given a complex ICA method:
%           >>    W = complexICA(data,options);
%           you can define a handle as:
%           >>    f = @(x) complexICA(x,options);
%           and feed the key-value pair ...,'ComplexICAMethod',f,...
%           into pwcica. Note, you may need to write a separate wrapper
%           function if complexICA has multiple outputs.
%     -> 'TimeInvariant' , ([1] , 0)
%           If set to 0, utulizes a transformation that scales the 
%           difference (tangent) vectors according to the 'SamplingRate' 
%           Default set to 1; in this case the mapping V is time-invariant
%           and is a complex linear structure.
%     -> 'SamplingRate', ([1], any real)
%           The sampling rate of observations. Set 'TimeInvariant' to 0 in
%           order to implement a sampling rate.
%     -> 'EpochLength',  ([K], an integer divisible by N if K = 1)
%           Define epoch length in case data should be epoched but is input
%           as stream. Ignored if K ~= 1. 
%     -> 'GPU' ([0] , 1)
%           Turn on to use MATLAB CUDA optimized 'ebm' complex ICA.
%     -> 'OptObject' (['pow3'] 'hyvar')
%           Specify optimization objective for complex 'fastICA'
%     -> 'Level' ([1] , any natural number > 0)
%           Width of time interval for sum/difference mapping into complex
%           phase space.
%     -> 'SGSmoothing' ([0] , 1)
%           Turn on to use Savitzky-Golay smoothing for mapping time series
%           data to complex phase space.
%     -> 'SGWindow' ([11], any odd natural number)
%           Specifiy window for Savitzky-Golay smoothing.
%     -> 'SGOrder' ([5], any natural number > 0)
%           Specify order of polynomial for Savitzky-Golay filter to fit.
%     -> 'ComplexICAPath' ([working directory\ComplexICA], string, or []);
%           Optionally specify the directory that external complex ICA
%           algorithms are stored in. If not set, pwcica assumes that
%           algorithms are located in a subdirectory caled 'ComplexICA' of
%           pwcica's working directory. Value may be set to a string of the
%           full path name of the relevant directory. If value is empty,
%           pwcica will NOT add anything to MATLAB's search path.
%
%   OUTPUTS:
%     W - Demixing matrix, [n,n] array.
%     metas (optional) - A MATLAB struct containing the following "meta"
%     data from the PWC-ICA problem:
%       metas.complexW - The complex ICA demixing matrix on whitened data.
%       metas.SphC - The complex whitening transformation used in pwcica.
%       metas.options - The cell array of key-value pairs input in
%       varargin.
%
%   EXAMPLE:
%       If data is an [n,N,K] array of n-dimensional data divided into 
%       epochs of length N, then:
%
%       >>[W,metas] = pwcica(data)
%       returns the time-invariant pwc-ica demixing matrix using a 
%       step-size of 1, derived with complex FastICA.
%
%       >>[W,metas] = pwcica(data,'Level',4);
%       returns the time-invariant pwc-ica demixing matrix using a
%       step-size of 4, derived with complex FastICA.
%
%       >>[W,metas] = pwcica(data,'SGSmoothing',1,'SGWindow',5,'SGOrder',3)
%       returns pwc-ica using an S-G filter of order 3, window size 5.
%
%       >>[W,metas] =
%       pwcica(data,'TimeInvariant',0,'SamplingRate',64,'ComplexICAMethod','ebm')
%       returns pwc-ica demixing matrix using step-size 1, a sampling rate
%       of 64 Hz, and Entropy Bound Minimization Complex ICA.

%   
%   **!! Important note on ComplexICA packages !!**
%
%   Complex ICA must be performed separately using one of the following
%   externally availabe functions:
%
%   Code to be used with pwcica:
%   FicaCPLX, by PETR TICHAVSKÝ
%   complex_ICA_EBM, by Xi-Lin Li
%   complex_ICA_EBM_GPU, by Xi-Lin Li, modified by Kenneth Ball
%   roubstica, by Vicente Zarzoso
%
%   Place the relevant function in an active MATLAB path so that pwcica can access the method. 
%
%
%   Author: Kenneth Ball, Postdoc; UTSA Dept of CS/ARL TNB-HRED
%   Copywrite: Kenneth Ball, 2015
%
%   Reference: Ball, K. R. et al. [2015] PWC-ICA: A Method for Ordered Blind 
%   Source Separation with Application to EEG. In Internal Review. 
%
%   If you find this code useful in your own analyses, please cite:
%
%   1. The PWC-ICA paper:
%   Ball, K. R., Bigdley-Shamlo, N., Mullins, T., Robbins, K.
%   PWC-ICA: A Method for Ordered Blind Source Separation with Application
%   to EEG. (To appear)
%
%   2. The paper associated with the Complex ICA method you use, especially
%   if you use one of the default packages. The default method 'FicaCPLX'
%   method's reference is:
%   
%   Koldovský, Z., & Tichavský, P. (2007). Blind Instantaneous Noisy 
%   Mixture Separation with Best Interference-Plus-Noise Rejection. 
%   In M. Davies, C. James, S. Abdallah, & M. Plumbley (Eds.), 
%   Independent Component Analysis and Signal Separation SE  - 91 
%   (Vol. 4666, pp. 730–737). Springer Berlin Heidelberg.
% 
% %
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA

function [ W , varargout ] = pwcica( data , varargin)

% Test to see if EEGLAB EEG struct is input as data:

if isa(data,'struct')
    EEG = data;
    data = EEG.data;
end

% Set relevant constants:
n = size(data,1);
N = size(data,2);

% Set Default Values:
if size(data,3) > 1
    % Data is input as epochs
    epochLength = N;
    epochCount = size(data,3);
    N = epochLength*epochCount;
else
    epochLength = N;
    epochCount = 1;
end
timeFlag = 0;
if exist('EEG','var')
    if isfield(EEG,'srate')
        dT = 1./EEG.srate;
%         timeFlag = 1;
    else
        dT = 1;
    end
else
    dT = 1;
end
complexICAMethod = 'fastICA';
gpu = 0;
level = 1;
% multiLevelFlag = 0;
% projectionType = 'average';
optimizationFlag = 'pow3';
savGolFlag = 0;
savGolWindow = 11;
savGolOrder = 5;
realFlag = 0;
WinitFlag = 0;
complexICAPath = mfilename('fullpath');
complexICAPath = [complexICAPath(1:end-7),'\ComplexICA'];

% Adjust for optional inputs
if (nargin> 1)
    if (nargin> 1 && rem(nargin,2) == 0)
        if length(varargin) == 1
            varargin = varargin{1};
        else
            fprintf('pwcica(): Optional key and value pairs do not match.')
            return
        end
    end
    
    for ii = 1:2:length(varargin)
        key = varargin{ii};
        val = varargin{ii+1};
        if strcmp(key,'EpochLength')
            if length(data,3) == 1
                epochLength = val;
                epochCount = N/val;
                if rem(N,val) ~= 0
                    fprintf('pwcica(): Total number of time frames must be divisible by Epoch Length.')
                    return
                end
            end
        elseif strcmp(key,'SamplingRate')
%             timeFlag = 1;
            if isempty(val)
                val = 1;
                timeFlag = 0;
            end
            dT = 1./val;
        elseif strcmp(key,'ComplexICAMethod')
            complexICAMethod = val;
        elseif strcmp(key,'GPU')
            gpu = val;
        elseif strcmp(key,'Level')
            level = val;
        elseif strcmp(key,'OptObject')
            optimizationFlag = val;
        elseif strcmp(key,'SGSmoothing')
            savGolFlag = val;
        elseif strcmp(key,'SGWindow')
            if rem(val,2)
                savGolWindow = val;
            else
                error('SG Window size must be odd');
            end
        elseif strcmp(key,'SGOrder')
            savGolOrder = val;
        elseif strcmp(key,'Real')
            realFlag = val;
        elseif strcmp(key,'Initializer')
            if val == 0
                WinitFlag = 0;
            else
                WinitFlag = 1;
                Winit = val;
            end
        elseif strcmp(key,'TimeInvariant')
            timeFlag = ~val;
        elseif strcmp(key,'ComplexICAPath')
            if isempth(val);
                complexICAPath = [];
            else
                complexICAPath = val;
            end
        end
    end
end


% Add paths for complex ICA algorithms:
if ~isempty(complexICAPath);
    addpath(genpath(complexICAPath));
end

% If input data is epoched, transform to stream first:
if size(data,3) > 1
    data = reshape(data,n,N);
end

% Zero-Mean the channels in the stream, in case data not already
% zero-meaned:
data1 = bsxfun(@minus,data,mean(data,2));
data = data1;

% Transform back to epochs:
data = reshape(data,n,epochLength,epochCount);

if savGolFlag == 0

    % Transform data epochs to pairwise form:
    
    data = cat(1,data(:,1:(end-level),:),data(:,(1+level):end,:));

    % Return to data stream form:

    data = reshape(data,2*n,(epochLength-level)*epochCount);

    % Generate tangent space transform and apply to pairwise data stream:
    if timeFlag
        V = [ .5*eye(n) , .5*eye(n) ; -1/(level*dT)*eye(n) , 1/(level*dT)*eye(n) ];
    else
        V = 1/sqrt(2)*[eye(n) , eye(n) ; -eye(n) , eye(n)];
    end

    data = V*data;
    
    % Zero Mean the data again.
    data = bsxfun(@minus,data,mean(data,2));

    % Transform into complex phase space:

    data = data(1:n,:) + 1i*data((n+1):(2*n),:);

    % Double the data in case it isn't already:

    data = double(data);
    
    % For real test case, take only the real-valued or base-space data:
    if realFlag
        data = real(data);
    end

    % Whiten in complex space. In this case, SphC is a transformation that
    % weakly whitens the data, in the sense that the covariance matrix of
    % the transformed data is the identity matrix, but Not the
    % pseudo-covariance matrix.
    
    SphC = inv(sqrtm((data*data')./size(data,2)));
    data = SphC*data;

    % Perform the complex ICA transformation:
    switch complexICAMethod
        case 'ebm'
            if gpu
                [complexW , ~ , ~] = complex_ICA_EBM_GPU(data);
            else
                [complexW , ~ , ~] = complex_ICA_EBM(data);
            end
        case 'robust'
            [~,~,~,complexW] = robustica(data,[],1e-6,1028,0,'o',0,[],1);
        case 'fastICA'
            if ~WinitFlag
                [complexW,~] = FicaCPLX(data,optimizationFlag);
            else
                [complexW,~] = FicaCPLX(data,optimizationFlag,true,Winit);
            end
        case 'infomax'
            [W,S] = runica(data,'Extended',1);
            complexW = W*S;
        otherwise
            complexW = complexICAMethod(data);
    end
    
else
    
    % Find the Savinsky-Golay filter for specified order and window size
    [~,g] = sgolay(savGolOrder,savGolWindow);
    halfWin = ((savGolWindow+1)/2)-1;
    complexData = [];
    if timeFlag
        for ii = (halfWin+1):(epochLength-halfWin)
            complexData = cat(2,complexData,...
                sum(bsxfun(@times,g(:,1)',data(:,(ii-halfWin):(ii+halfWin),:)),2) ...
                + 1i*sum(bsxfun(@times,g(:,2)',data(:,(ii-halfWin):(ii+halfWin),:)),2)/dT ...
                );
        end
    else
        for ii = (halfWin+1):(epochLength-halfWin)
        complexData = cat(2,complexData,...
            sum(bsxfun(@times,g(:,1)',data(:,(ii-halfWin):(ii+halfWin),:)),2) ...
            + 1i*sum(bsxfun(@times,g(:,2)',data(:,(ii-halfWin):(ii+halfWin),:)),2) ...
            );
        end
    end
    
    data = double(complexData);
    clear complexData;
    
    data = reshape(data,n,(epochLength-halfWin*2)*epochCount);
    
    % Zero Mean the data again.
    data = bsxfun(@minus,data,mean(data,2));
    
    % For real test case, take only the real-valued or base-space data:
    if realFlag
        data = real(data);
    end
    
    % Whiten in complex space. In this case, SphC is a transformation that
    % weakly whitens the data, in the sense that the covariance matrix of
    % the transformed data is the identity matrix, but Not the
    % pseudo-covariance matrix.
    
    SphC = inv(sqrtm((data*data')./size(data,2)));
    data = SphC*data;

    % Perform the complex ICA transformation:
    switch complexICAMethod
        case 'ebm'
            if gpu
                [complexW , ~ , ~] = complex_ICA_EBM_GPU(data);
            else
                [complexW , ~ , ~] = complex_ICA_EBM(data);
            end
        case 'robust'
            [~,~,~,complexW] = robustica(data,[],1e-6,1028,0,'o',0,[],1);
        case 'fastICA'
            if ~WinitFlag
                [complexW,~] = FicaCPLX(data,optimizationFlag);
            else
                [complexW,~] = FicaCPLX(data,optimizationFlag,true,Winit);
            end
        case 'infomax'
            [W,S] = runica(data,'Extended',1);
            complexW = W*S;
        otherwise
            complexW = complexICAMethod(data);
    end
    
end

W = real(complexW*SphC);

scaling = repmat(sqrt(mean(inv(W).^2))',[1 size(W,1)]);
% feats = (scaling.*W)*data1;
variance = var(((scaling.*W)*data1).');
[~,varIndex] = sort(variance,'descend');
W = W(varIndex,:);

metas.complexW = complexW;
metas.SphC = SphC;
metas.options = varargin;
metas.varIndex = varIndex;

varargout{1} = metas;



end