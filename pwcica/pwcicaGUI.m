% pwcicaGUI -- Calls GUI for pwcica within EEGLAB
%
% Usage: options = pwcicaGUI(EEG)
% Meant to be called within EEGLAB GUI environment, from the
% eegplugin_pwcica environment.
% Arguments:
% INPUT: Accepts an EEG type struct from eeglab.
% OUTPUT: options, a cell array of key-value pairs to be fed into
% pop_pwcica(), and subsequently pwcica().
%
%   Author: Kenneth Ball, Postdoc; UTSA Dept of CS/ARL TNB-HRED
%   Copywrite Kenneth Ball, 2015
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

function options = pwcicaGUI(EEG)


if isfield(EEG,'srate')
%     srate = num2str(EEG.srate);
    srate = EEG.srate;
else
    srate = num2str('');
end

thegeometry = { [.5 .5 .5] [1 .5] [1 .5] [1 .5] [.2 .15 1.15] [.2 .15 1.15] [.2 .15 1.15] [1 .5] };
theuilist = { ...
    {'style', 'text', 'string' , '_____________'},...
    {'style', 'text', 'FontWeight', 'bold', 'string', 'PWC-ICA'},...
    {'style', 'text', 'string' , '_____________'},...
    {'style', 'text', 'string', 'Time invariant?'},...
    {'style', 'checkbox', 'value', 1, 'tag', 'ti'},...
    {'style', 'text', 'string', 'Step-size:    '},...
    {'style', 'edit', 'string', 1, 'tag','level'},...
    {'style', 'text', 'string', 'Complex ICA Method:'},...
    {'style', 'edit', 'string', 'fastICA', 'tag', 'method'},...
    {'style', 'text', 'string', ' '},...
    {'style', 'text', 'string', 'Ex:'},...
    {'style', 'text', 'string', '[fastICA , ebm , robustica]'},...
    {'style', 'text', 'string', ' '},...
    {'style', 'text', 'string', 'Or:'},...
    {'style', 'text', 'string', 'Function handle of the form:'},...
    {'style', 'text', 'string', ' '},...
    {'style', 'text', 'string', ' '},...
    {'style', 'text', 'string', '>> W = fun(ComplexData);'},...
    {'style', 'text', 'string', 'SamplingRate:   '},...
    {'style', 'text', 'string', [num2str(srate),' Hz']}...
    };

[~ , ~ , ~ , result] = inputgui( 'geometry', thegeometry, 'uilist', theuilist, ...
    'helpcom', 'pophelp(''pop_pwcica'');', 'title', 'Run PWC-ICA -- pop_pwcica()');


% result = inputgui('geometry' , { 1 [1 .2] [1 .2] [1 .5] [1 .2] },...
%     'uilist', { ...
%     {'Style', 'text', 'string', 'PWC-ICA'},...
%     {'Style', 'text', 'string', 'Time invariant?'},...
%     {'Style', 'checkbox', 'string', ' ', 'value', 1, 'tag', 'ti'},...
%     {'Style', 'text', 'string', 'Step-size = '},...
%     {'Style', 'edit', 'string', '1', 'tag','level'},...
%     {'Style', 'text', 'Complex ICA Method:'},...
%     {'Style', 'edit', 'string', 'fastica', 'tag', 'method'},...
%     {'Style', 'text', 'string', 'SamplingRate = '},...
%     {'Style', 'edit', 'string', srate, 'tag', 'srate'}...
%     });

options = [];
if isempty(result)
    return
end
if isfield(result,'ti')
    if ~isempty(result.ti)
        options = cat(2,options,{'TimeInvariant',result.ti});
    end
end
if isfield(result,'level')
    if ~isempty(result.level)
        if isnumeric(result.level)
            options = cat(2,options,{'Level',result.level});
        elseif ischar(result.level)
            options = cat(2,options,{'Level',str2num(result.level)});
        end
    end
end
if isfield(result,'method')
    if ~isempty(result.method) && ischar(result.method)
        options = cat(2,options,{'ComplexICAMethod',result.method});
    end
end

% Check to see if sampling rate is unspecified for time variant pwcica
if isempty(srate) && ~result.ti 
    thegeometry1 = { 1 1 [.6 .3 .1] };
    theuilist1 = { ...
        {'style', 'text', 'string', 'Warning: No sampling rate specified in EEG struct, but time-variant pwcica selected.'},...
        {'style', 'text', 'string', 'Please specify sampling rate'},...
        {'style', 'text', 'string', 'Sampling Rate     ='},...
        {'style', 'edit', 'string', '', 'tag', 'srate'},...
        {'style', 'text', 'string', ' Hz'}
    };
    [~ , ~ , ~ , result1] = inputgui( 'geometry', thegeometry1, 'uilist', theuilist1, ...
        'helpcom', 'pophelp(''pop_pwcica'');', 'title', 'Warning: Specify Sampling Rate');
    
    if isfield(result1,'srate')
        if ~isempty(result1.srate)
            if isnumeric(result1.srate)
                options = cat(2,options,{'SamplingRate',result1.srate});
            elseif ischar(result1.srate)
                options = cat(2,options,{'SamplingRate',str2num(result1.srate)});
            end
        else
            display('Error: You must set a sampling rate to perform time-variant pwcica.');
            options = [];
            return
        end
    end
    
    if isempty(result1)
        options = [];
        return
    end
    
elseif ~isempty(srate)
    options = cat(2,options,{'SamplingRate',srate});
end





end