% demo to show the usage of complex_ICA_EBM
clear
close all;

        % QAM sources
        T = 500;           % the sample size
        N = 10;             % generate 10 QAM sources
        for n = 1 : N
            M = 2^n;        % M is the order of QAM, here from 2 to 1024
            x = [0:M-1];    % x the real part symbols
            y = qammod(x,M);    % y is the final symbols, real and imaginary parts
            s(n,:) = randsrc(1,T,[y; ones(1,M)/M]); % generate QAM sources, each symbols have the same probability
            if M == 2 
                noise = awgn(s(n,:),20,'measured') - s(n,:);    % for real BPSK, added noise is real
                noise = noise.*exp(2*pi*sqrt(-1)*rand(size(noise)));     % random rotate to get complex noise
                s(n,:) = s(n,:) + noise;                        % add complex noise;
            else
                s(n,:) = awgn(s(n,:),20,'measured');    % by default, add complex Gaussian noise;
            end
 
            % scatterplot(s(n,:))     % show the scatter plot;
            
        end


% remove mean and power normalization
for n = 1 : N
    s(n,:) = s(n,:) - mean(s(n,:));
    s(n,:) = sqrt(T)*s(n,:)/norm(s(n,:));
end
        
% mixing
A = randn(N,N)+sqrt(-1)*randn(N,N);
x = A*s;

% separation
W = complex_ICA_EBM(x);

% showing results
T = abs(W*A);
for n=1:N
    T(n,:)=T(n,:)/(max(abs(T(n,:))));
end
T=kron(T,ones(50,50));
figure;
imshow(abs(T))



