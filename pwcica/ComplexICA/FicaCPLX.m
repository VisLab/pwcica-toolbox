function [W1U,Wsym]=FicaCPLX(X, g, SaddleTest, ini);
% 
%  complex FastICA separation algorithm with a seddle point test.
%  The algorithm proceeds by (1) performing a symmetric complex FastICA,
%  (2) saddle point test to reduce probability of side local minima, and
%  (3) 1 unit refinement for each separated component.
%
%  Inputs:
%  =======
%  X.....input signal dxN, where N is number of samples
%  g.....used nonlinearity: pow3 ... recommended for finite alphabet signals
%                           hyvar ... proposed by Hyvarinen
%  SaddleTest ... the test can be turned off
%  ini .......... initialization of the algorithm
%
%  Outputs
%  =======
%  W1U .... estimated de-mixing matrix produced by the one unit refinement
%  Wsym ... estimated demixing matrix produced by the symmetric FastICA.
%
%  NOTE: The separated components are defined here as S_1U = W1U*X
%                                                     S_sym = Wsym*X
%        unlike the convention in Bingham & Hyvarinen !!! 
%
[dim N]=size(X);
X=X-mean(X,2)*ones(1,N);  %  removing sample mean
jj=sqrt(-1);
%Default values of parameters
if nargin<4
    ini=eye(dim);
end
if nargin<3
    SaddleTest=true;
end
if nargin<2
    g='pow3';
end
epsilon=0.00001; %Stop criterion
fineepsilon=1e-5; %Stop criterion for post-estimation
repeat=1;
rot2d=[1/sqrt(2) 1/sqrt(2);-1/sqrt(2) 1/sqrt(2)];
MaxIt=100; %% Maximum number of FastICA iterations
MaxItAfterSaddleTest=30; %%Maximum number of iterations after a saddle point was indicated
FicaMaxIt=50; %%Maximum number of improving iterations
min_correlation=0.8; %additive noise...0.75, noise free... 0.95, turn off (unstable)...0
%
C = cov(X');
CC = C^(-1/2);
Z = CC*X;
%CC=eye(2);
%FastICA
W=ini;
W_before_decor=W;
W=W*real(inv(W'*W)^(1/2));
Wold=zeros(dim,dim);
crit=zeros(1,dim);
NumIt=0;
%%% Symmetric approach
while repeat
while (1-min(crit)>epsilon && NumIt<MaxIt)% && sum(double(changed>10))<2)
   Wold=W;
   switch g
    case 'pow3'
%      W=(Z*((Z'*W).^2.*conj(Z'*W)))/N-2*W;
     u = Z'*W; W=(Z*((u).^2.*conj(u)))/N-2*W;
    case 'hyvar'
    % W=(Z*((Z'*W).^2.*conj(Z'*W)))/N-2*W;
     Z1=W'*Z;
     g1=1./(1+abs(Z).^2);
     p1=Z*(conj(Z1.*g1)')/N;
     p2=mean(g1-(abs(Z1).^2).*(g1.^2),2);
     W=p1-W.*(ones(dim,1)*conj(p2'));  
   end
  % W_before_decor=W;
   W=W*(W'*W)^(-1/2);
   crit=sum(abs(W.*Wold));
   NumIt=NumIt+1; 
end %while iteration
fprintf('FicaCPLX: NumIt= %d\n',NumIt);
repeat=0;
%%%
%%%      doing dim of 2-dimensional real-valued FastICA to
%%%      separate real and imaginary parts in each component
%%%
for i=1:dim;
    W(:,i)=rotuj(W(:,i),Z);
end      
%%%%  NOW, TEST OF SADDLE POINTS
if SaddleTest
  SaddleTest=false; %%The SaddleTest may be done only one times
  u=Z'*W;
  table1=(mean((abs(u).^2).^2)-2).^2;
  rotated=zeros(1,dim);
  checked=1:dim;
  for i=checked
    for j=checked(checked>i)
      if (~rotated(i) && ~rotated(j))
        h=[u(:,i) u(:,j)]*rot2d;
        ctrl=sum((mean((abs(h).^2).^2)-2).^2); 
        if ctrl>table1(i)+table1(j)
            %bad extrem indicated
            rotated(i)=1;rotated(j)=1; %do not test the rotated signals anymore
            repeat=1;
            W(:,[i j])=W(:,[i j])*rot2d;  
            W(:,i)=rotuj(W(:,i),Z);
            W(:,j)=rotuj(W(:,j),Z);
            crit=zeros(1,dim);
            repeat=1;
        end
    end  % of if rotated
    end  % of j checked
  end    % of i checked
end      % of SaddleTest
end      % of repeat
Wsym=W;
%%%                 NOW, CONTINUE TO COMPUTE 1FICA      
for i=1:dim
 w=Wsym(:,i);
 w=w/norm(w);
 winit=w;
 wold=zeros(dim,1);
 NumIt=0;
 while (abs(w'*wold)<1-epsilon & NumIt<FicaMaxIt & abs(w'*winit)>min_correlation)
    wold=w;
    z=w'*Z;
    switch g
      case 'pow3'
        g2=abs(z).^2;
        p1=Z*((conj(z).*g2).')/N;
        p2=2;             % p2=2*mean(g2);     
      case 'hyvar'
        g2=1./(1+abs(z).^2);
        p1=Z*((conj(z).*g2).')/N;
        p2=mean(g2+(abs(z).^2).*(-g2.^2)); 
    end 
    w=p1-w*p2;     
    w=w/norm(w);
    NumIt=NumIt+1;
 end
 if abs(w'*winit)>min_correlation
    W1U(:,i)=w;
 else
    W1U(:,i)=winit;
 end  
end
Wsym=(CC*Wsym)';
W1U=(CC*W1U)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function rotuj - real FastICA for 2 signals
function w=rotuj(w,Z)
%
% finds such a complex multiplication factor such that the real and imaginary
% parts of w'*Z are as independent as possible (here the independence is measured by the kurtosis
% criterion)
%
  [dim,N]=size(Z);
  jj=sqrt(-1);
  epsilon=0.0001; %Stop criterion
  repeat=1;
  rot2d=[1/sqrt(2) 1/sqrt(2);-1/sqrt(2) 1/sqrt(2)];
  MaxIt=100; %% Maximum number of FastICA iterations
  xcompl=w'*Z;
  X2=[real(xcompl); imag(xcompl)];
  C2=X2*X2'/N;
  CC2=C2^(-1/2);
  Z2=CC2*X2;
  W2=eye(2);
  crit2=zeros(1,dim); NumIt2=0;
  while repeat  %% performing two-dimensional symmetric FastICA
    while (1-min(crit2)>epsilon && NumIt2<MaxIt)
          Wold2=W2;
          W2=(Z2*((Z2'*W2).^ 3))/N-3*W2;
          W2=W2*(W2'*W2)^(-1/2);
          crit2=sum(abs((W2.*Wold2)));
          NumIt2=NumIt2+1;
    end
    repeat=0;
    u=Z2'*W2;  %%% now, do the test of saddle points
    table1=(mean((u.^2).^2)-3).^2;
    u=u*rot2d;
    table2=(mean((u.^2).^2)-3).^2;
    if sum(table2)>sum(table1)
       W2=W2*rot2d;
       repeat=1;
    end
  end  % of repeat    
  Wsymm2=W2*CC2;
  RW1=real(w)*Wsymm2(1,1)+imag(w)*Wsymm2(1,2);
  RW2=real(w)*Wsymm2(2,1)+imag(w)*Wsymm2(2,2);
  w=RW1+jj*RW2;
  w=w/norm(w);