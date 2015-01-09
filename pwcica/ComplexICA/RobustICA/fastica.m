function [S, H, iter, W] = fastica(X, kurtsign, tol, max_it, prewhi, deftype, dimred, Wini, verbose)
  
% Kurtosis-based FastICA method for deflationary ICA/BSS of:
%
% A. Hyvarinen and E. Oja, "A fast fixed-point algorithm for independent component analysis", 
% Neural Computation, vol. 9, no. 7, pp. 1483-149, Oct. 1997.
%
% A. Hyvarinen, "Fast and robust fixed-point algorithms for independent component analysis",
% IEEE Transactions on Neural Networks, vol. 10, no. 3, pp. 626-634, May 1999.
%
%
% SYNTAX: [S, H, iter, W] = fastica(X, kurtsign, tol, max_it, prewhi, deftype, dimred, Wini, verbose);
%
%
% OUTPUTS:
%         S       : estimated sources signals (one row per signal)
%
%         H       : estimated mixing matrix
%
%         iter    : number of iterations (one element per extracted source)
%
%         W       : final extracting vectors found in FastICA iteration.
%
%
% INPUTS:
%
%         X       : observed signals (one row per signal)
%
%         kurtsign: not used (included for compatibility with RobustICA)
%
%         tol     : threshold for statistically-significant termination test equivalent to
%                        ||wn - w||/||w|| < tol/sqrt(sample size);
%                   termination test not used if tol < 0
%
%         max_it  : maximum number of iterations per extracting vector
%
%         prewhi  : prewhitening (via SVD of the observed data matrix);
%                   no prewhitening if prewhi = 0
%
%         deftype : deflation type: 'o'rthogonalization, 'r'egression
%
%         dimred  : dimensionality reduction in regression if parameter different from zero;
%                   not used (included for compatibility with RobustICA)
%
%         Wini    : extracting vector initialization for FastICA iterative search
%
%         verbose : verbose operation if different from 0.        
%
%
% HISTORY:
%
% -- 2010/02/16: Version 2 release ----------------------------------------------------------------
%
% 2010/02/10: - created and added to RobustICA package by Vicente Zarzoso
%               (University of Nice Sophia Antipolis, France).


[n, T] = size(X);


%%% remove mean

X = X - mean(X')'*ones(1, T);


%%%%% Prewhitening (if required)

if prewhi
    if verbose,  disp(' '); disp('>>> Prewhitening'); end
    [V, D, U] = svd(X', 0);     % economy SVD of data matrix
    B = sqrt(T)*U*diag(1./diag(D));    % whitening matrix, such that Z = B'*X
    Z = sqrt(T)*V';             % PCA source estimate
else
    B = eye(n);
    Z = X;
end;


%%%%% FastICA algorithm

W = zeros(n);   % extracting vectors
I = eye(n);     
P = I;          % projection matrix for deflationary orthogonalization (if required)

if ~exist('Wini'), Wini = I; % initialization
else 
    if isempty(Wini), Wini = I; end; 
end;

tol = tol/sqrt(T);          %  a statistically-significant termination threshold 
tol2 = sign(tol)*tol^2/2;   % obtain an equivalent termination test for FastICA iteration
iter = zeros(1, n);         % number of iterations                                        

if deftype == 'r'
    do_reg = 1;
    def_text = 'regression-based deflation';
elseif deftype == 'o'
    do_reg = 0;
    def_text = 'deflationary orthogonalization'; 
else
    error('Deflation method not recognized (''deftype'' parameter): use ''r'' or ''o''.')
end

if verbose, disp(['>>> FastICA with ', def_text]), end


for k = 1:n

if verbose, disp(['> src#', num2str(k)]), end

it = 0;
keep_going = 1;

w = P*Wini(:, k);   % initialization
w = w/norm(w);      % normalization


while keep_going

    it = it + 1;        % increase number of iterations
    
    y = w'*Z;
    wn = P*(w - (Z*(y.^3)'/(3*T)));     % equivalent to wn = E[z(w'z)^3] - 3w ...                            
    wn = wn/norm(wn);                   % ... if followed by normalization
   
    %check for convergence
    th = abs(1 - abs(wn'*w));

    w = wn;
   
    if th < tol2 || it >= max_it || any(isnan(wn))    
            % finish when convergence or too many iterations or divergence       
            keep_going = 0;    
    end 
    
end %while keep_going

W(:, k) = w;    % estimated extracting vector
s = w'*Z;       % estimated source
S(k, :) = s;
iter(k) = it;   % number of  iterations

if verbose, disp([num2str(it), ' iterations']), end

if do_reg
    Z = deflation_regression(Z, s); % regression + subtraction
else
    P = I - W*W';   % projection matrix for orthogonalization (if required)
end

end % for k


if verbose, disp(['Total number of iterations: ', num2str(sum(iter))]), end


%%%%% estimate mixing matrix

if prewhi
    H = B*W;      
else
    H = inv(W');
end
