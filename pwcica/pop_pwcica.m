% pop_pwcica() - Calls pwcica for EEGLAB EEG data structures.
%
% Usage:
%   >>  OUTEEG = pop_pwcica( INEEG, type, param3 );
%
% Inputs:
%   INEEG   - input EEG dataset
%   type    - type of processing. 1 process the raw
%             data and 0 the ICA components.
%   param3  - additional parameter: a [1,2n] cell array of n key-value pairs.
%    
% Outputs:
%   OUTEEG  - output dataset
%
% See also:
%   SAMPLE, EEGLAB 

% Copyright (C) <year>  <name of author>
%
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
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [EEG, com] = pop_pwcica( EEG, typeproc, param3 )

% the command output is a hidden output that does not have to
% be described in the header

com = ''; % this initialization ensure that the function will return something
          % if the user press the cancel button            

% display help if not enough arguments
% ------------------------------------
if nargin < 2
	help pop_pwcica;
	return;
end;	

% pop up window
% -------------
if nargin < 3
	promptstr    = { 'Enter the parameter:' };
	inistr       = { '0' };
	result       = inputdlg( promptstr, 'Title of window', 1,  inistr);
	if length( result ) == 0 return; end;

	param3   	 = eval( [ '[' result{1} ']' ] ); % the brackets allow to process matlab arrays
end;

% call function sample either on raw data or ICA data
% ---------------------------------------------------
if typeproc == 1
    addSRateFlag = 0;
    if isfield(EEG,'srate')
        addSRateFlag = 1;
    end
	for ii = 1:2:length(param3)
        if strcmp(param3{ii},'SamplingRate')
            addSRateFlag = 0;
            break;
        end
    end
    if addSRateFlag
        param3 = cat(2,param3,{'SamplingRate',EEG.srate});
    end
    W = pwcica(EEG.data,param3);
    scaling = repmat(sqrt(mean(inv(W).^2))',[1 size(W,1)]);
    display('Scaling components to RMS microvolt');
    W = scaling.*W;
    EEG.icaweights = W;
    EEG.icasphere = eye(size(EEG.data,1));
    EEG.icawinv = inv(W);
else
    error('Wrong procedure type');
	if ~isempty( EEG.icadata )
		pwcica( EEG.icadata );
	else
		error('You must run ICA first');
	end;	
end;	 

% return the string command
% -------------------------
foo = '{ ';
if ischar(param3{2})
    foo = [foo,'''',param3{1}, ''', ''',param3{2},''''];
else
    foo = [foo,'''',param3{1}, ''', ',num2str(param3{2})];
end
for ii = 3:2:length(param3);
    if ischar(param3{ii+1})
        foo = [foo, ', ''', param3{ii}, ''', ''', param3{ii+1},''''];
    else
        foo = [foo, ', ''', param3{ii}, ''', ', num2str(param3{ii+1})];
    end
end
foo = [foo,' }'];
    
% com = sprintf('pop_pwcica( %s, %d, [%s] );', inputname(1), typeproc, int2str(param3));
com = sprintf('pop_pwcica( %s, %d, %s );', inputname(1), typeproc, foo);

return;
