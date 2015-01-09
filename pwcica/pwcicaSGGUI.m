function options = pwcicaSGGUI(EEG)

if isfield(EEG,'srate')
    srate = EEG.srate;
else
    srate = num2str('');
end

thegeometry = { 1 [1 .2] [1 .5] 1 [1 .5] [1 .5] [1 .5]};
theuilist = { ...
    {'style', 'text', 'string', 'PWC-ICA w/ Savitzky-Golay Filtering'},...
    {'style', 'text', 'string', 'Time invariant?'},...
    {'style', 'checkbox', 'value', 1, 'tag', 'ti'},...
    {'style', 'text', 'string', 'Complex ICA Method:'},...
    {'style', 'edit', 'string', 'fastICA', 'tag', 'method'},...
    {'style', 'text', 'string', 'Options: [fastICA , ebm , robustica]'},...
    {'style', 'text', 'string', 'SamplingRate     = '},...
    {'style', 'text', 'string', [num2str(srate),' Hz']}...
    {'style', 'text', 'string', 'SG Filter Order  = '},...
    {'style', 'edit', 'string', 3, 'tag','order'},...
    {'style', 'text', 'string', 'SG Window Length = '},...
    {'style', 'edit', 'string', 7, 'tag','window'},...
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

options = {'SGSmoothing',1};
if isempty(result)
    return
end

if isfield(result,'ti')
    if ~isempty(result.ti)
        options = cat(2,options,{'TimeInvariant',result.ti});
    end
end
if isfield(result,'method')
    if ~isempty(result.method) && ischar(result.method)
        options = cat(2,options,{'ComplexICAMethod',result.method});
    end
end
% if isfield(result,'srate')
%     if ~isempty(result.srate)
%         if isnumeric(result.srate)
%             options = cat(2,options,{'SamplingRate',result.srate});
%         elseif ischar(result.srate)
%             options = cat(2,options,{'SamplingRate',str2num(result.srate)});
%         end
%     end
% end
if isfield(result,'order')
    if ~isempty(result.order)
        if isnumeric(result.order)
            options = cat(2,options,{'SGOrder',result.order});
        elseif ischar(result.order)
            options = cat(2,options,{'SGOrder',str2num(result.order)});
        end
    else
        display('Warning: No Savitzky-Golay Order specified, default will be used');
    end
end
if isfield(result,'window')
    if ~isempty(result.window)
        if isnumeric(result.window)
            options = cat(2,options,{'SGWindow',result.window});
        elseif ischar(result.order)
            options = cat(2,options,{'SGWindow',str2num(result.window)});
        end
    else
        display('Warning: No Savitzky-Golay Window length specified, default will be used');
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