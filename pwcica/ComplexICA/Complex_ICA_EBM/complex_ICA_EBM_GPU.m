function [W, Ahat, Shat] = complex_ICA_EBM_GPU( X )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% complex ICA-EBM: complex ICA by entropy bound minimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input: mixtures X
% Outputs: demixing matrix W, estimations of mixing matrix and sources (Ahat, Shat)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reference
% Xi-Lin Li and Tulay Adali, "Complex independent component analysis by
% entropy bound minimization," to appear in IEEE Transaction on Circuit and Systems I. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program by Xi-Lin Li
% Please contact me at lixilin@umbc.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tolerance = 1e-4;   % stopping condition 
show_cost = 0;

[N,T] = size(X);
[Xc, P] = pre_processing( X ); Xc = gpuArray(Xc);

% load 8 nonlinearities, but we only use 4 of them
K_real = 8;
load nf_table
real_nf1 = nf1;     % x^4
real_nf2 = nf2;     % |x|
real_nf3 = nf3;     % |x|/(1+|x|)
real_nf4 = nf4;     % |x|/(1+x^2)
real_nf5 = nf5;     % x|x|/(10+|x|)
real_nf6 = nf6;     % x/(1+|x|)
real_nf7 = nf7;     % x/(1+x^2)
real_nf8 = nf8;

% here begins the sea algorithm to provide the initial guess
W = randn(N)+sqrt(-1)*randn(N);
W = inv(sqrtm(W*W'))*W;  % W = gpuArray(W);
last_W = W;
C = Xc*Xc.'/T;      % C is the pseudo covariance matrix
maxiter_sea = 100;
Cost = zeros(maxiter_sea,1);
min_cost = inf;
cost_increase_counter = 0;
max_cost_increase_number = 10;
for iter = 1 : maxiter_sea
    
    for n = 1 : N
        v = gpuArray(W(n,:).');
        y = v.'*Xc;
        %v = conj(Xc)*( y.*y.*conj(y) ).'/T - 2*v - (v.'*C*v)*conj(C*v);
        foo = -v.';
        foo = foo*C;
        foo = foo*v;
        bar = C*v;
        bar = conj(bar);
        foo = -2*v + foo*bar;
        bar = conj(Xc)*( y.*y.*conj(y) ).'/T ;
        bar = gather(bar);
        v = bar + foo;
        
        W(n,:) = gather(v.');
        
        % evaluate the cost
        z = y;
        z_real = real(z);
        z_imag = imag(z);
        sigma_R2 = gather(sum( z_real.^2 )/T);
        sigma_I2 = gather(sum( z_imag.^2 )/T);
        sigma_R = sqrt(sigma_R2);
        sigma_I = sqrt(sigma_I2);
        rho = gather(sum( z_real.*z_imag )/T);
        Delta1 = sigma_R2*sigma_I2 - rho^2;
        u = z_real/sigma_R;
        v = sigma_R*z_imag/sqrt(Delta1) - rho*z_real/sigma_R/sqrt(Delta1);
        uu = u.*u;
        vv = v.*v;
        sign_u = sign(u);
        abs_u = sign_u.*u;
        sign_v = sign(v);
        abs_v = sign_v.*v;
        
        Cost(iter) = Cost(iter) + 0.5*log(Delta1);
        
        NE_Boundu = zeros(K_real,1);
        NE_Boundv = zeros(K_real,1);
        EGu = zeros(K_real,1);
        EGv = zeros(K_real,1);
        
        % G1 = x^4
        EGu(1) = gather(sum(uu.*uu))/T; ;%EGu(1) = sum(uu.*uu)/T;
        if EGu(1)<real_nf1.min_EGx
            NE_Boundu(1) = simplified_ppval( real_nf1.pp, real_nf1.min_EGx );
        else if EGu(1)>real_nf1.max_EGx
                NE_Boundu(1) =  simplified_ppval( real_nf1.pp, real_nf1.max_EGx );
            else
                NE_Boundu(1) =  simplified_ppval( real_nf1.pp, EGu(1) );
            end
        end
       
        % G3 = |x|/(1+|x|)
        EGu(3) = gather(sum( abs_u./(1+abs_u) ))/T;
        if EGu(3)<real_nf3.min_EGx
            NE_Boundu(3) = simplified_ppval( real_nf3.pp_slope, real_nf3.min_EGx ) * ( EGu(3) - real_nf3.min_EGx );
            NE_Boundu(3) = simplified_ppval( real_nf3.pp, real_nf3.min_EGx ) + abs( NE_Boundu(3) ); 
        else if EGu(3)>real_nf3.max_EGx
                NE_Boundu(3) =  simplified_ppval( real_nf3.pp_slope, real_nf3.max_EGx ) * ( EGu(3) - real_nf3.max_EGx );
                NE_Boundu(3) =  simplified_ppval( real_nf3.pp, real_nf3.max_EGx ) + abs( NE_Boundu(3) );
            else
                NE_Boundu(3) =  simplified_ppval( real_nf3.pp, EGu(3) );
            end
        end
        
        % G5 = x*|x|/(10+|x|)
        EGu(5) = gather(sum( u.*abs_u./(10+abs_u) ))/T;
        if EGu(5)<real_nf5.min_EGx
            NE_Boundu(5) = simplified_ppval( real_nf5.pp_slope, real_nf5.min_EGx ) * ( EGu(5) - real_nf5.min_EGx );
            NE_Boundu(5) = simplified_ppval( real_nf5.pp, real_nf5.min_EGx ) + abs( NE_Boundu(5) ); 
        else if EGu(5)>real_nf5.max_EGx
                NE_Boundu(5) =  simplified_ppval( real_nf5.pp_slope, real_nf5.max_EGx ) * ( EGu(5) - real_nf5.max_EGx );
                NE_Boundu(5) =  simplified_ppval( real_nf5.pp, real_nf5.max_EGx ) + abs( NE_Boundu(5) );
            else
                NE_Boundu(5) =  simplified_ppval( real_nf5.pp, EGu(5) );
            end
        end
        
        % G7 = x/(1+x^2)
        EGu(7) = gather(sum( u./(1+u.*u) ))/T;
        if EGu(7)<real_nf7.min_EGx
            NE_Boundu(7) = simplified_ppval( real_nf7.pp_slope, real_nf7.min_EGx ) * ( EGu(7) - real_nf7.min_EGx );
            NE_Boundu(7) = simplified_ppval( real_nf7.pp, real_nf7.min_EGx ) + abs( NE_Boundu(7) ); 
        else if EGu(7)>real_nf7.max_EGx
                NE_Boundu(7) =  simplified_ppval( real_nf7.pp_slope, real_nf7.max_EGx ) * ( EGu(7) - real_nf7.max_EGx );
                NE_Boundu(7) =  simplified_ppval( real_nf7.pp, real_nf7.max_EGx ) + abs( NE_Boundu(7) );
            else
                NE_Boundu(7) =  simplified_ppval( real_nf7.pp, EGu(7) );
            end
        end
               
        % G1 = x^4
        Gv(1) = gather(sum(vv.*vv))/T; % EGv(1) = sum(vv.*vv)/T;
        if EGv(1)<real_nf1.min_EGx
            NE_Boundv(1) = simplified_ppval( real_nf1.pp, real_nf1.min_EGx );
        else if EGv(1)>real_nf1.max_EGx
                NE_Boundv(1) =  simplified_ppval( real_nf1.pp, real_nf1.max_EGx );
            else
                NE_Boundv(1) =  simplified_ppval( real_nf1.pp, EGv(1) );
            end
        end
       
        % G3 = |x|/(1+|x|)
        EGv(3) = gather(sum( abs_v./(1+abs_v) ))/T;
        if EGv(3)<real_nf3.min_EGx
            NE_Boundv(3) = simplified_ppval( real_nf3.pp_slope, real_nf3.min_EGx ) * ( EGv(3) - real_nf3.min_EGx );
            NE_Boundv(3) = simplified_ppval( real_nf3.pp, real_nf3.min_EGx ) + abs( NE_Boundv(3) );
        else if EGv(3)>real_nf3.max_EGx
                NE_Boundv(3) =  simplified_ppval( real_nf3.pp_slope, real_nf3.max_EGx ) * ( EGv(3) - real_nf3.max_EGx );
                NE_Boundv(3) =  simplified_ppval( real_nf3.pp, real_nf3.max_EGx ) + abs( NE_Boundv(3) );
            else
                NE_Boundv(3) =  simplified_ppval( real_nf3.pp, EGv(3) );
            end
        end
        
        % G5 = x|x|/(10+|x|)
        EGv(5) = gather(sum( v.*abs_v./(10+abs_v) ))/T;
        if EGv(5)<real_nf5.min_EGx
            NE_Boundv(5) = simplified_ppval( real_nf5.pp_slope, real_nf5.min_EGx ) * ( EGv(5) - real_nf5.min_EGx );
            NE_Boundv(5) = simplified_ppval( real_nf5.pp, real_nf5.min_EGx ) + abs( NE_Boundv(5) );
        else if EGv(5)>real_nf5.max_EGx
                NE_Boundv(5) =  simplified_ppval( real_nf5.pp_slope, real_nf5.max_EGx ) * ( EGv(5) - real_nf5.max_EGx );
                NE_Boundv(5) =  simplified_ppval( real_nf5.pp, real_nf5.max_EGx ) + abs( NE_Boundv(5) );
            else
                NE_Boundv(5) =  simplified_ppval( real_nf5.pp, EGv(5) );
            end
        end
        
        % G7 = x/(1+x^2)
        EGv(7) = gather(sum( v./(1+v.*v) ))/T;
        if EGv(7)<real_nf7.min_EGx
            NE_Boundv(7) = simplified_ppval( real_nf7.pp_slope, real_nf7.min_EGx ) * ( EGv(7) - real_nf7.min_EGx );
            NE_Boundv(7) = simplified_ppval( real_nf7.pp, real_nf7.min_EGx ) + abs( NE_Boundv(7) );
        else if EGv(7)>real_nf7.max_EGx
                NE_Boundv(7) =  simplified_ppval( real_nf7.pp_slope, real_nf7.max_EGx ) * ( EGv(7) - real_nf7.max_EGx );
                NE_Boundv(7) =  simplified_ppval( real_nf7.pp, real_nf7.max_EGx ) + abs( NE_Boundv(7) );
            else
                NE_Boundv(7) =  simplified_ppval( real_nf7.pp, EGv(7) );
            end
        end
  
        Cost(iter) = Cost(iter) - max(NE_Boundu) - max(NE_Boundv);
        
    end
    
    if Cost(iter) < min_cost
        min_cost = Cost(iter);
        cost_increase_counter = 0;
        best_W = last_W;
    else
        cost_increase_counter = cost_increase_counter + 1;
    end

    W = inv(sqrtm(W*W'))*W;
    if cost_increase_counter > max_cost_increase_number
        break;
    end
    if 1-min(abs(diag(W*last_W'))) < tolerance
        break;
    else
        last_W = W;
    end
end
W = best_W;
if show_cost
    figure;
    subplot(1,3,1)
    plot(Cost(1:iter))
    title('Initialization')
end

max_cost_increase_number = 5;
%first use stochastic gradient search algorithm
% W = complex_ICA_EBM_north( Xc, W, 1000, 1/5, max_cost_increase_number, 0, 1, show_cost );
W = complex_ICA_EBM_north( Xc, W, 1000, gpuArray(1/5), max_cost_increase_number, 0, 1, show_cost );
%then refine the solution
% W = complex_ICA_EBM_north( Xc, W, 500, 1/25, 3, 1, 0, show_cost );       % refinement
W = complex_ICA_EBM_north( Xc, W, 500, gpuArray(1/25), 3, 1, 0, show_cost );       % refinement
W = gather(W);

W = W*P;
% [EOF]

% if you want to get the sources, it is
Ahat = inv(W);
Shat = W*X;















function W = complex_ICA_EBM_north( X, W0, max_iter_north, mu0_north, max_cost_increase_number, type2, stochastic_search, show_cost )
% complex ICA-EBM north version
% assume the mixtures have been normalized
% if type2~=0, the second kind of entropy bound will be used too
% if stochastic_search~=0, then do the stochastic gradient search
[N,T] = size(X);
R_xxt = X*X.'/T;

% load the nonlinearities. but we only 4 of them
K_real = 8;
load nf_table
real_nf1 = nf1;     % x^4
real_nf2 = nf2;     % |x|
real_nf3 = nf3;     % |x|/(1+|x|)
real_nf4 = nf4;     % |x|/(1+x^2)
real_nf5 = nf5;     % x|x|/(10+|x|)
real_nf6 = nf6;     % x/(1+|x|)
real_nf7 = nf7;
real_nf8 = nf8;

% we only use 2 of them
if type2
    K_complex = 4;
    load complex_nf_table
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Nonorthogonal ICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error = 0;
cost_increase_counter = 0;
mu = mu0_north;
W = W0;
best_W = W;
last_W = W;
Cost = zeros(max_iter_north,1);
min_cost = inf;
min_cost_queue = min_cost*ones(max_iter_north,1);
max_negentropy = zeros(N,1);
negentropy_array = zeros(N,1);
for iter = 1 : max_iter_north
   
    Cost(iter) = gather(-log(abs(det(W'*W))));
    %W = W(randperm(N),:);
    
    for n = 1 : N
        
        if N>7
            
            if n==1
                Wn = W(2:N,:);
                inv_Q = inv( Wn*Wn' );
            else
                n_last = n-1;
                Wn_last = [ W(1:n_last-1,:);W(n_last+1:N,:) ];
                w_current = W(n,:)';
                w_last = W(n_last,:)';
                c = Wn_last*( w_last - w_current );
                c(n_last) = 0.5*( w_last'*w_last - w_current'*w_current );
                e_last = zeros(N-1,1);
                e_last(n_last) = 1;
            
                temp1 = inv_Q*c;
                temp2 = inv_Q(:,n_last);
                inv_Q_plus = inv_Q - temp1*temp2'/(1+temp1(n_last));
            
                temp1 = inv_Q_plus'*c;
                temp2 = inv_Q_plus(:,n_last);
                inv_Q = inv_Q_plus - temp2*temp1'/(1+c'*temp2);
                % inv_Q is Hermitian
                inv_Q = (inv_Q+inv_Q')/2;
            end

            temp1 = randn(N, 1);
            W_n = [ W(1:n-1,:);W(n+1:N,:) ];
            h = temp1 - W_n'*inv_Q*W_n*temp1;
            
        else
            temp1 = randn(N, 1);
            temp2 = [W(1:n-1,:); W(n+1:N,:)];
            h = temp1 - temp2'*inv(temp2*temp2')*temp2*temp1;
        end
        
        w = gpuArray(W(n,:)');
        z = w'*X;
        z_real = real(z);
        z_imag = imag(z);
        sigma_R2 = gather(sum( z_real.^2 ))/T;
        sigma_I2 = gather(sum( z_imag.^2 ))/T;
        sigma_R = sqrt(sigma_R2);
        sigma_I = sqrt(sigma_I2);
        rho = gather(sum( z_real.*z_imag ))/T;
        Delta1 = sigma_R2*sigma_I2 - rho^2;
        u = z_real/sigma_R;
        v = sigma_R*z_imag/sqrt(Delta1) - rho*z_real/sigma_R/sqrt(Delta1);
        uu = u.*u;
        vv = v.*v;
        sign_u = sign(u);
        abs_u = sign_u.*u;
        sign_v = sign(v);
        abs_v = sign_v.*v;
        
        if type2
            Delta2 = sigma_I2*z_real.^2 + sigma_R2*z_imag.^2 - 2*rho*z_real.*z_imag;
            rr = Delta2/Delta1;
            r = sqrt(rr);
        end
        
        Cost(iter) = Cost(iter) + 0.5*log(gather(Delta1));
        negentropy_array(n) = -0.5*log(gather(Delta1));         % we consider the negentropy
        
        NE_Boundu = zeros(K_real,1);
        NE_Boundv = zeros(K_real,1);
        EGu = zeros(K_real,1);
        EGv = zeros(K_real,1);
        
        % G1 = x^4
        EGu(1) = gather(sum(uu.*uu))/T;
        if EGu(1)<real_nf1.min_EGx
            NE_Boundu(1) = simplified_ppval( real_nf1.pp, real_nf1.min_EGx );
        else if EGu(1)>real_nf1.max_EGx
                NE_Boundu(1) =  simplified_ppval( real_nf1.pp, real_nf1.max_EGx );
            else
                NE_Boundu(1) =  simplified_ppval( real_nf1.pp, EGu(1) );
            end
        end
       
        % G3 = |x|/(1+|x|)
        EGu(3) = gather(sum( abs_u./(1+abs_u) ))/T;
        if EGu(3)<real_nf3.min_EGx
            NE_Boundu(3) = simplified_ppval( real_nf3.pp_slope, real_nf3.min_EGx ) * ( EGu(3) - real_nf3.min_EGx );
            NE_Boundu(3) = simplified_ppval( real_nf3.pp, real_nf3.min_EGx ) + abs( NE_Boundu(3) ); 
        else if EGu(3)>real_nf3.max_EGx
                NE_Boundu(3) =  simplified_ppval( real_nf3.pp_slope, real_nf3.max_EGx ) * ( EGu(3) - real_nf3.max_EGx );
                NE_Boundu(3) =  simplified_ppval( real_nf3.pp, real_nf3.max_EGx ) + abs( NE_Boundu(3) );
            else
                NE_Boundu(3) =  simplified_ppval( real_nf3.pp, EGu(3) );
            end
        end
        
        % G5 = x*|x|/(10+|x|)
        EGu(5) = gather(sum( u.*abs_u./(10+abs_u) ))/T;
        if EGu(5)<real_nf5.min_EGx
            NE_Boundu(5) = simplified_ppval( real_nf5.pp_slope, real_nf5.min_EGx ) * ( EGu(5) - real_nf5.min_EGx );
            NE_Boundu(5) = simplified_ppval( real_nf5.pp, real_nf5.min_EGx ) + abs( NE_Boundu(5) ); 
        else if EGu(5)>real_nf5.max_EGx
                NE_Boundu(5) =  simplified_ppval( real_nf5.pp_slope, real_nf5.max_EGx ) * ( EGu(5) - real_nf5.max_EGx );
                NE_Boundu(5) =  simplified_ppval( real_nf5.pp, real_nf5.max_EGx ) + abs( NE_Boundu(5) );
            else
                NE_Boundu(5) =  simplified_ppval( real_nf5.pp, EGu(5) );
            end
        end
        
        % G7 = x/(1+x^2)
        EGu(7) = gather(sum( u./(1+u.*u) ))/T;
        if EGu(7)<real_nf7.min_EGx
            NE_Boundu(7) = simplified_ppval( real_nf7.pp_slope, real_nf7.min_EGx ) * ( EGu(7) - real_nf7.min_EGx );
            NE_Boundu(7) = simplified_ppval( real_nf7.pp, real_nf7.min_EGx ) + abs( NE_Boundu(7) ); 
        else if EGu(7)>real_nf7.max_EGx
                NE_Boundu(7) =  simplified_ppval( real_nf7.pp_slope, real_nf7.max_EGx ) * ( EGu(7) - real_nf7.max_EGx );
                NE_Boundu(7) =  simplified_ppval( real_nf7.pp, real_nf7.max_EGx ) + abs( NE_Boundu(7) );
            else
                NE_Boundu(7) =  simplified_ppval( real_nf7.pp, EGu(7) );
            end
        end
               
        % G1 = x^4
        EGv(1) = gather(sum(vv.*vv))/T;
        if EGv(1)<real_nf1.min_EGx
            NE_Boundv(1) = simplified_ppval( real_nf1.pp, real_nf1.min_EGx );
        else if EGv(1)>real_nf1.max_EGx
                NE_Boundv(1) =  simplified_ppval( real_nf1.pp, real_nf1.max_EGx );
            else
                NE_Boundv(1) =  simplified_ppval( real_nf1.pp, EGv(1) );
            end
        end
       
        % G3 = |x|/(1+|x|)
        EGv(3) = gather(sum( abs_v./(1+abs_v) ))/T;
        if EGv(3)<real_nf3.min_EGx
            NE_Boundv(3) = simplified_ppval( real_nf3.pp_slope, real_nf3.min_EGx ) * ( EGv(3) - real_nf3.min_EGx );
            NE_Boundv(3) = simplified_ppval( real_nf3.pp, real_nf3.min_EGx ) + abs( NE_Boundv(3) );
        else if EGv(3)>real_nf3.max_EGx
                NE_Boundv(3) =  simplified_ppval( real_nf3.pp_slope, real_nf3.max_EGx ) * ( EGv(3) - real_nf3.max_EGx );
                NE_Boundv(3) =  simplified_ppval( real_nf3.pp, real_nf3.max_EGx ) + abs( NE_Boundv(3) );
            else
                NE_Boundv(3) =  simplified_ppval( real_nf3.pp, EGv(3) );
            end
        end
        
        % G5 = x|x|/(10+|x|)
        EGv(5) = gather(sum( v.*abs_v./(10+abs_v) ))/T;
        if EGv(5)<real_nf5.min_EGx
            NE_Boundv(5) = simplified_ppval( real_nf5.pp_slope, real_nf5.min_EGx ) * ( EGv(5) - real_nf5.min_EGx );
            NE_Boundv(5) = simplified_ppval( real_nf5.pp, real_nf5.min_EGx ) + abs( NE_Boundv(5) );
        else if EGv(5)>real_nf5.max_EGx
                NE_Boundv(5) =  simplified_ppval( real_nf5.pp_slope, real_nf5.max_EGx ) * ( EGv(5) - real_nf5.max_EGx );
                NE_Boundv(5) =  simplified_ppval( real_nf5.pp, real_nf5.max_EGx ) + abs( NE_Boundv(5) );
            else
                NE_Boundv(5) =  simplified_ppval( real_nf5.pp, EGv(5) );
            end
        end
        
        % G7 = x/(1+x^2)
        EGv(7) = gather(sum( v./(1+v.*v) ))/T;
        if EGv(7)<real_nf7.min_EGx
            NE_Boundv(7) = simplified_ppval( real_nf7.pp_slope, real_nf7.min_EGx ) * ( EGv(7) - real_nf7.min_EGx );
            NE_Boundv(7) = simplified_ppval( real_nf7.pp, real_nf7.min_EGx ) + abs( NE_Boundv(7) );
        else if EGv(7)>real_nf7.max_EGx
                NE_Boundv(7) =  simplified_ppval( real_nf7.pp_slope, real_nf7.max_EGx ) * ( EGv(7) - real_nf7.max_EGx );
                NE_Boundv(7) =  simplified_ppval( real_nf7.pp, real_nf7.max_EGx ) + abs( NE_Boundv(7) );
            else
                NE_Boundv(7) =  simplified_ppval( real_nf7.pp, EGv(7) );
            end
        end
             
        p = 0;
        q = 0;
        p0 = 0;
        q0 = 0;
        bound_real = -inf;
        for p = 1 : K_real
            for q = 1 : K_real
                if NE_Boundu(p) + NE_Boundv(q) > bound_real
                    bound_real = NE_Boundu(p) + NE_Boundv(q);
                    p0 = p;
                    q0 = q;
                end
            end
        end
        
        if type2
            
            NE_Boundr = zeros(K_complex,1);
            EGr = zeros(K_complex,1);

            % G1 = r^4
            EGr(1) = gather(sum(rr.*rr))/T;
            if EGr(1)<nf1.min_EGr
                NE_Boundr(1) =  simplified_ppval( nf1.pp, nf1.min_EGr );
            else if EGr(1)>nf1.max_EGr
                    NE_Boundr(1) =  simplified_ppval( nf1.pp, nf1.max_EGr );
                else
                    NE_Boundr(1) =  simplified_ppval( nf1.pp, EGr(1) );
                end
            end

            % G3 = r/(1+r)
            EGr(3) = gather(sum( r./(1+r) ))/T;
            if EGr(3)<nf3.min_EGr
                NE_Boundr(3) =  simplified_ppval( nf3.pp_slope, nf3.min_EGr ) * ( EGr(3) - nf3.min_EGr );
                NE_Boundr(3) =  simplified_ppval( nf3.pp, nf3.min_EGr ) + abs( NE_Boundr(3) );
            else if EGr(3)>nf3.max_EGr
                    NE_Boundr(3) =  simplified_ppval( nf3.pp_slope, nf3.max_EGr ) * ( EGr(3) - nf3.max_EGr );
                    NE_Boundr(3) =  simplified_ppval( nf3.pp, nf3.max_EGr ) + abs( NE_Boundr(3) );
                else
                    NE_Boundr(3) =  simplified_ppval( nf3.pp, EGr(3) );
                end
            end

            [max_NEr, max_i] = max( NE_Boundr );  
            
        end
        
        if type2 == 0 || bound_real > max_NEr
            
            Cost(iter) = Cost(iter) - bound_real;
            negentropy_array(n) = negentropy_array(n) + bound_real;
        
            grad_delta1 = ( 0.5*(sigma_I2 - sigma_R2) + sqrt(-1)*rho )*R_xxt*conj(w);
            
            if stochastic_search
                weight = gpuArray.rand(1,T);
            else
                weight = gpuArray.ones(1,T);
            end
            sumWeight = sum(weight);

            switch p0

                case 1
                    % G = x^4
                    % g = 4x^3
                    vEGu = simplified_ppval( real_nf1.pp_slope, EGu(1) );
                    grad = 0.5*grad_delta1/Delta1;
                    grad = grad - vEGu*X*( 4*weight.*uu.*u )'/sumWeight/2/sigma_R;
                    grad = grad - vEGu*sum( -4*weight.*uu.*u.*z_real )/sumWeight/4/sigma_R^3*R_xxt*conj(w);

                case 3
                    % G = |x|/(1+|x|)
                    % g = sign(x)/(1+|x|)^2
                    vEGu = simplified_ppval( real_nf3.pp_slope, min(max(EGu(3), real_nf3.min_EGx), real_nf3.max_EGx) );
                    grad = 0.5*grad_delta1/Delta1;
                    grad = grad - vEGu*X*( weight.*sign_u./(1+abs_u).^2 )'/sumWeight/2/sigma_R;
                    grad = grad - vEGu*sum( -weight.*sign_u./(1+abs_u).^2.*z_real )/sumWeight/4/sigma_R^3*R_xxt*conj(w);
                    
                case 5
                    % G = x|x|/(10+|x|)
                    % g = |x|(20+|x|)/(10+|x|)^2
                    vEGu = simplified_ppval( real_nf5.pp_slope, min(max(EGu(5), real_nf5.min_EGx), real_nf5.max_EGx) );
                    grad = 0.5*grad_delta1/Delta1;
                    grad = grad - vEGu*X*( weight.*abs_u.*(20+abs_u)./(10+abs_u).^2 )'/sumWeight/2/sigma_R;
                    grad = grad - vEGu*sum( -weight.*abs_u.*(20+abs_u)./(10+abs_u).^2.*z_real )/sumWeight/4/sigma_R^3*R_xxt*conj(w);
                    
                case 7
                    % G = x/(1+x^2)
                    % g = (1-x^2)/(1+x^2)^2
                    vEGu = simplified_ppval( real_nf7.pp_slope, min(max(EGu(7), real_nf7.min_EGx), real_nf7.max_EGx) );
                    grad = 0.5*grad_delta1/Delta1;
                    grad = grad - vEGu*X*( weight.*(1-u.*u)./(1+u.*u).^2 )'/sumWeight/2/sigma_R;
                    grad = grad - vEGu*sum( -weight.*(1-u.*u)./(1+u.*u).^2.*z_real )/sumWeight/4/sigma_R^3*R_xxt*conj(w);

            end
            
            if stochastic_search
                weight = gpuArray.rand(1,T);%weight = rand(1,T);
            else
                weight = gpuArray.ones(1,T);%weight = ones(1,T);
            end

            switch q0

                case 1
                    % G = x^4
                    vEGv = simplified_ppval( real_nf1.pp_slope, EGv(1) );
                    grad = grad + vEGv*X*( 4*weight.*vv.*v )'/sumWeight/2/sigma_R/sqrt(Delta1)*(rho+sqrt(-1)*sigma_R2);
                    grad = grad + vEGv*sum( 4*weight.*vv.*vv )/sumWeight/2/Delta1*grad_delta1;
                    grad = grad - vEGv*sum( 4*weight.*vv.*v.*( z_imag + (rho/sigma_R2+2*sqrt(-1))*z_real ) )/sumWeight/4/sigma_R/sqrt(Delta1)*R_xxt*conj(w);

                case 3
                    % G = |x|/(1+|x|)
                    % g = sign(x)/(1+|x|)^2
                    vEGv = simplified_ppval( real_nf3.pp_slope, min(max(EGv(3), real_nf3.min_EGx), real_nf3.max_EGx) );
                    grad = grad + vEGv*X*( weight.*sign_v./(1+abs_v).^2 )'/sumWeight/2/sigma_R/sqrt(Delta1)*(rho+sqrt(-1)*sigma_R2);
                    grad = grad + vEGv*sum( weight.*abs_v./(1+abs_v).^2 )/sumWeight/2/Delta1*grad_delta1;
                    grad = grad - vEGv*sum( weight.*sign_v./(1+abs_v).^2.*( z_imag + (rho/sigma_R2+2*sqrt(-1))*z_real ) )/sumWeight/4/sigma_R/sqrt(Delta1)*R_xxt*conj(w);

                case 5
                    % G = x|x|/(10+|x|)
                    % g = |x|(20+|x|)/(10+|x|)^2
                    vEGv = simplified_ppval( real_nf5.pp_slope, min(max(EGv(5), real_nf5.min_EGx), real_nf5.max_EGx) );
                    grad = grad + vEGv*X*( weight.*abs_v.*(20+abs_v)./(10+abs_v).^2 )'/sumWeight/2/sigma_R/sqrt(Delta1)*(rho+sqrt(-1)*sigma_R2);
                    grad = grad + vEGv*sum( weight.*abs_v.*(20+abs_v)./(10+abs_v).^2.*v )/sumWeight/2/Delta1*grad_delta1;
                    grad = grad - vEGv*sum( weight.*abs_v.*(20+abs_v)./(10+abs_v).^2.*( z_imag + (rho/sigma_R2+2*sqrt(-1))*z_real ) )/sumWeight/4/sigma_R/sqrt(Delta1)*R_xxt*conj(w);
                    
                case 7
                    % G = x/(1+x^2)
                    % g = (1-x^2)/(1+x^2)^2
                    vEGv = simplified_ppval( real_nf7.pp_slope, min(max(EGv(7), real_nf7.min_EGx), real_nf7.max_EGx) );
                    grad = grad + vEGv*X*( weight.*(1-v.*v)./(1+v.*v).^2 )'/sumWeight/2/sigma_R/sqrt(Delta1)*(rho+sqrt(-1)*sigma_R2);
                    grad = grad + vEGv*sum( weight.*(1-v.*v)./(1+v.*v).^2.*v )/sumWeight/2/Delta1*grad_delta1;
                    grad = grad - vEGv*sum( weight.*(1-v.*v)./(1+v.*v).^2.*( z_imag + (rho/sigma_R2+2*sqrt(-1))*z_real ) )/sumWeight/4/sigma_R/sqrt(Delta1)*R_xxt*conj(w);
                    
            end
            
        else
            
            Cost(iter) = Cost(iter) - max_NEr;
            negentropy_array(n) = negentropy_array(n) + max_NEr;

            grad_delta1 = ( 0.5*(sigma_I2 - sigma_R2) + sqrt(-1)*rho )*R_xxt*conj(w);
            grad_delta2_part1 = sigma_I2*z_real - sqrt(-1)*sigma_R2*z_imag + sqrt(-1)*rho*z ;
            grad_delta2_part2 = 0.5*(z_imag.^2 - z_real.^2) + sqrt(-1)*z_real.*z_imag;
            inv_r = 1./(r+eps);
            
            if stochastic_search
                weight = gpuArray.rand(1,T);%weight = rand(1,T);
            else
                weight = gpuArray.ones(1,T);%weight = ones(1,T);
            end
            sumWeight = sum(weight);
            
            switch max_i
            
                case 1
                    % G = r^4
                    % g = 4*r^3
                    vEGr = simplified_ppval( nf1.pp_slope, EGr(1) );
                    grad = 0.5*(  1/Delta1 + 1/Delta1^2*vEGr*sum(4*weight.*rr.*Delta2)/sumWeight  )*grad_delta1;
                    grad = grad - 0.5*1/Delta1*vEGr* X*( 4*weight.*rr.*grad_delta2_part1 ).'/sumWeight;
                    grad = grad  - 0.5*1/Delta1*vEGr* sum( 4*weight.*rr.*grad_delta2_part2 )/sumWeight*R_xxt*conj(w);
                
                case 3
                    % G = r/(1+r)
                    vEGr = simplified_ppval( nf3.pp_slope, min(max(EGr(3), nf3.min_EGr), nf3.max_EGr) );
                    grad = 0.5*(  1/Delta1 + 1/Delta1^2*vEGr*sum(weight.*inv_r./(1+r).^2.*Delta2)/sumWeight  )*grad_delta1;
                    grad = grad - 0.5*1/Delta1*vEGr* X*( weight.*inv_r./(1+r).^2.*grad_delta2_part1 ).'/sumWeight;
                    grad = grad  - 0.5*1/Delta1*vEGr* sum( weight.*inv_r./(1+r).^2.*grad_delta2_part2 )/sumWeight*R_xxt*conj(w);
                                
            end

        end

        grad = grad - h/(w'*h);
        grad = grad - (real(w'*grad))*w;
        grad = grad / norm(grad);
        
        w1 = w - mu*grad;
        w1 = w1 / norm(w1);
        W(n,:) = gather(w1');
        
    end
    
    if Cost(iter) < min_cost
        min_cost = Cost(iter);
        best_W = last_W;
        max_negentropy = negentropy_array;
        cost_increase_counter = 0;
    else
        cost_increase_counter = cost_increase_counter + 1;
    end
    min_cost_queue(iter) = min_cost;
    
    if cost_increase_counter > max_cost_increase_number
        if stochastic_search
            if mu>1/20
                mu = mu/2;
                cost_increase_counter = 0;
                W = best_W;
                last_W = W;
                continue;
            else
                break;
            end
        else
            if mu>1/200
                mu = mu/2;
                cost_increase_counter = 0;
                W = best_W;
                last_W = W;
                continue;
            else
                break;
            end
        end
    end

    last_W = W;
    
end
W = best_W;

% sort all the component
[max_negentropy, index_sort] = sort(max_negentropy, 'descend');
W = W(index_sort, :);

if show_cost
    if stochastic_search
        subplot(1,3,2);
    else
        subplot(1,3,3);
    end
    hold on; plot([1:iter], Cost([1:iter]));
    %hold on; plot([1:iter], min_cost_queue([1:iter]),'r-');
    xlabel('Number of iterations')
    ylabel('Cost')
    title('Nonorthogonal ICA')
end
%[EOF]



function [X, P] = pre_processing( X )
% pre-processing program
[N,T] = size(X);
% remove DC
Xmean=mean(X,2);
X = X - Xmean*ones(1,T);    

% spatio pre-whitening 1 
R = X*X'/T;                 
P = inv_sqrtmH(R);  %P = inv(sqrtm(R));
X = P*X;



function A = inv_sqrtmH( B )
% 
[V,D] = eig(B);
d = diag(D);
d = 1./sqrt(d);
A = V*diag(d)*V';


function v = simplified_ppval(pp,xs)
% a simplified version of ppval

b = pp.breaks;
c = pp.coefs;
l = pp.pieces;
k = 4;  % k = pp.order;
dd = 1; % dd = pp.dim;

% find index
index=0;
middle_index=0;
if xs>b(l)
    index = l;
else if xs<b(2)
        index = 1;
    else 
        
        low_index = 1;
        high_index = l;
        while 1
            
            middle_index = round( (low_index + high_index)/2 );
            if b( middle_index ) > xs
                high_index = middle_index;
            else
                low_index = middle_index;
            end
            
            if low_index == high_index-1
                index = low_index;
                break;
            end
            
        end
        
    end
end
    

% now go to local coordinates ...
xs = xs-b(index);

% ... and apply nested multiplication:
   v = c(index,1);
   for i=2:k
      v = xs.*v + c(index,i);
   end