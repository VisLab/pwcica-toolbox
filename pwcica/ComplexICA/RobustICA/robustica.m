function [S, H, iter, W] = robustica(X, kurtsign, tol, max_it, prewhi, deftype, dimred, Wini, verbose)
 
% Kurtosis-based RobustICA method for deflationary ICA/BSS (see references below for details).
%
%
% SYNTAX: [S, H, iter, W] = robustica(X, kurtsign, tol, max_it, prewhi, deftype, dimred, Wini, verbose);
%
%
% OUTPUTS:
%         S       : estimated sources signals (one row per signal, one column per sample)
%
%         H       : estimated mixing matrix
%
%         iter    : number of iterations (one element per extracted source)
%
%         W       : estimated extracting vectors
%                  (acting on whitened observations if prewhitened is required;
%                   otherwise, acting on given observations).
%
%
% INPUTS:
%         X       : observed signals (one row per signal, one column per sample)
%
%         kurtsign: source kurtosis signs (one element per source or empty);
%                   maximize absolute normalized kurtosis if element = 0;
%                   if empty, maximize absolute normalized kurtosis for all sources
%
%         tol     : threshold for statistically-significant termination test of the type
%                         ||wn - p*w||/||w|| < tol/sqrt(sample size);   (up to a phase shift p)
%                   termination is also tested by comparing the gradient norm according to: 
%                         ||g|| < tol/sqrt(sample size);
%                   termination test is not used if tol < 0, so that the algorithm runs the maximum
%                   number of iterations (except if optimal step size is null)
%
%         max_it  : maximum number of iterations per extracted source
%
%         prewhi  : prewhitening (via SVD of the observed data matrix);
%                   no prewhitening if prewhi = 0
%
%         deftype : deflation type: 'o'rthogonalization, 'r'egression
%
%         dimred  : dimensionality reduction in regression if parameter different from zero;
%                   (not used in deflationary orthogonalization)
%
%         Wini    : extracting vectors initialization for RobustICA iterative search;
%                   if empty or not specified, identity matrix of suitable dimensions is used
%
%         verbose : verbose operation if parameter different from 0 (default: quiet operation).
%      
%
% EXAMPLES:
%
% - RobustICA with prewhitening, regression-based deflation, identity matrix initialization (default), 
%   up to 1000 iterations per source, termination threshold 1e-3/(sample size), without aiming at
%   any specific source (default):
%
%      >> S = robustica(X, [], 1e-3, 1e3, 1, 'r', 0);
%
% - RobustICA without prewhitening, with regression-based deflation:
%
%      >> S = robustica(X, [], 1e-3, 1e3, 0, 'r', 0);
%
% - RobustICA without prewhitening, with regression-based deflation, random initialization:
%
%      >> S = robustica(X, [], 1e-3, 1e3, 0, 'r', 0, randn(size(X, 1)));
%
% - RobustICA without prewhitening, with regression-based deflation and dimensionality reduction:
%
%      >> S = robustica(X, [], 1e-3, 1e3, 0, 'r', 1);
%
% - RobustICA with prewhitening, regression-based deflation, identity matrix initialization,
%   verbose operation:
%
%      >> S = robustica(X, [], 1e-3, 1e3, 1, 'r', 0, [], 1);
%
% - RobustICA with prewhitening, deflationary orthogonalization:
%
%      >> S = robustica(X, [], 1e-3, 1e3, 1, 'o', 0);
%
% - RobustICA with prewhitening, deflationary orthogonalization, and 10 iterations per independent
%   component:
%
%      >> S = robustica(X, [], -1, 10, 1, 'o', 0);
%
% - RobustICA with prewhitening, regression-based deflation, targeting first the sub-Gaussian
%   and then the super-Gaussian sources in a square mixture of 5 sub-Gaussian and 5 super-Gaussian sources:
%  
%      >> S = robustica(X, [ones(1,5), -ones(1,5)], 1e-3, 1e3, 1, 'r', 0);
%
% - RobustICA with prewhitening, regression-based deflation, targeting first a sub-Gaussian source:
%  
%      >> S = robustica(X, [-1, zeros(1,size(X, 1)-1)], 1e-3, 1e3, 1, 'r', 0);
%
%
% REFERENCES:
%
% - V. Zarzoso and P. Comon, <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/tnn10.pdf">"Robust independent component analysis by iterative maximization</a>
%   <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/tnn10.pdf">of the kurtosis contrast with algebraic optimal step size"</a>, 
%   IEEE Transactions on Neural Networks, vol. 21, no. 2, pp. 248-261, Feb. 2010.
%
% - V. Zarzoso and P. Comon, <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/ica07.pdf">"Comparative Speed Analysis of FastICA"</a>, 
%   in: Proceedings ICA-2007, 7th International Conference on Independent Component Analysis
%       and Signal Separation, London, UK, September 9-12, 2007, pp. 293-300.
% 
% - V. Zarzoso, P. Comon and M. Kallel,  <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/eusipco06.pdf">"How Fast is FastICA?"</a>, 
%   in: Proceedings EUSIPCO-2006, XIV European Signal Processing Conference, 
%       Florence, Italy, September 4-8, 2006. 
%
%
% Please, report any bugs, comments or suggestions to <a href = "mailto:zarzoso@i3s.unice.fr">zarzoso(a)i3s.unice.fr</a>.
%
%
% HISTORY:
% 
%    <Please add modification date here>: - <please add modification details here>
%       
% -- 2010/02/16: Version 2 release ----------------------------------------------------------------
%
%    2010/02/09: - added termination test based on gradient norm
%                
%    2009/03/02: - project extracting vector before normalization
% 
%    2009/02/02: - variable 'thmu' (for step-size based termination test) removed, as it was not used
%
% -- 2008/03/31: Version 1 release ----------------------------------------------------------------
%
%    2008/12/03: - modified help info about output parameter W
%                 (extracting vectors act on whitened observation if prewhitening is required)
% 
%    2008/03/26: - added this help
%
%    2008/03/13: - created by Vicente Zarzoso (University of Nice Sophia Antipolis, France).



[n, T] = size(X);


%%% remove mean

X = X - mean(X')'*ones(1, T);


if ~exist('verbose'), verbose = 0; end

%%%%% Prewhitening (if required)

if prewhi
    if verbose, disp(' '); disp('>>> Prewhitening'); end
    [V, D, U] = svd(X', 0);     % economy SVD of data matrix
    B = U*D/sqrt(T);            % PCA mixing-matrix estimate
    Z = sqrt(T)*V';             % PCA source estimate
else
    Z = X;
end;


%%%%% RobustICA algorithm

dimobs = n;     % number of remaining observations (may change under dimensionality reduction)

W = zeros(n);   % extracting vectors
I = eye(n);  
P = I;          % projection matrix for deflationary orthogonalization (if required)

if ~exist('Wini'), Wini = I; % initialization
else 
    if isempty(Wini), Wini = I; end; 
end;

tol = tol/sqrt(T);          % a statistically-significant termination threshold 
tol2 = sign(tol)*tol^2/2;   % the same threshold in terms of extracting vectors' absolute scalar product
iter = zeros(1, n);         % number of iterations                                        

if isempty(kurtsign), kurtsign = zeros(1, n); end
    
if deftype == 'r'
    do_reg = 1;
    typ_text = 'regression-based deflation';
    if dimred 
       typ_text = [typ_text, ' and dimensionality reduction'];     % only used in regression mode
   else
       typ_text = [typ_text, ', no dimensionality reduction'];   % default
    end
elseif deftype == 'o'
    do_reg = 0;
    typ_text = 'deflationary orthogonalization'; 
else
    error('Deflation method not recognized (''deftype'' parameter): use ''r'' or ''o''. Please type ''help robustica'' for details.')
end

if verbose, disp(' '); disp(['>>> RobustICA with ', typ_text]); end

    
% iterate over all sources
for k = 1:n

if verbose, disp(['> source #', num2str(k), ':']); end

it = 0;
keep_going = 1;

w = Wini(:, k);     % initialization
if do_reg, w = w((n-dimobs+1):n); end; % keep only required number of components

w = w/norm(w);      % normalization
w = P*w;            % project onto extracted vectors' orthogonal subspace (if deflationary orthogonalization)
                   
signkurt = kurtsign(k); % kurtosis sign of next source to be estimated


% iterate to extract one source
while keep_going

    it = it + 1;
    
    % compute KM optimal step size for gradient descent
    [g, mu_opt, norm_g] = kurt_gradient_optstep(w, Z, signkurt, P);  
         
    % update extracting vector
    wn = P*(w + mu_opt*g);      % update extracting vector and project if required 

    wn = wn/norm(wn);   % normalize 
   
    % extracting vector convergence test
    th = abs(1 - abs(wn'*w));
        
    w = wn;
    
    if th < tol2 || norm_g < tol || it >= max_it || mu_opt == 0 
        % finish when extracting vector converges, the gradient is too small, 
        % too many iterations have been run, or the optimal step-size is zero
       keep_going = 0;   
   end % if th < tol2 ...
     
end % while keep_going

if do_reg
    W(:,k) = [w; zeros(n-dimobs, 1)];  % zero-padding to account for dimensionality reduction in regression
else
    W(:, k) = w;    % estimated extracting vector
end % if do_reg

s = w'*Z;       % estimated source
S(k, :) = s;
iter(k) = it;   % number of  iterations

if verbose, disp([num2str(it), ' iterations']); end

if do_reg
    Z = deflation_regression(Z, s, dimred); % regression + subtraction
    dimobs = size(Z, 1);    % recompute observation dimension, just in case it has changed during deflation
    P = eye(dimobs);        % P is not required, but its dimensions should decrease according to Z
else
    P = I - W*W';     % projection matrix for orthogonalization (if required)
end % if do_reg

end % for k


if verbose, disp(' '); disp(['Total number of iterations: ', num2str(sum(iter))]); disp(' '); end


%%%%% Mixing matrix estimation

H = X*S'*pinv(S*S');  % LS estimate
