%%% Decoder Stochastic %%%
clear all;
% Parameter
EbN0 = 0:1:7;
EbN0_lin = 10.^(EbN0/10);
R = 4/7;
L = 7; VNj = 7; CNi = 3;
maxNumErrs = 180;
maxNum = 1e6;

%SD parameters
Nmax = 1;
LBernoulli = 1024;
BPiterMax = 20;

encTxF = [1 1 1 1 1 1 1]; %[all zero vector: BPSK] 

for i = 1:1:length(EbN0)
  totErr = 0;
  num = 0;

  while((totErr < maxNumErrs) && (num < maxNum)) 
    %Noise Addition

    sigma = 1/sqrt(2*EbN0_lin(i)*R);
    cgnoise = sigma.*randn(1,length(encTxF));
    %cgnoise = 0;
    y = encTxF + cgnoise;
    
    % Initialize (1)
    initVNj = zeros(VNj,LBernoulli);
    %AWGN
    for j = 1:VNj
      %channelInput = 2*y(j) ./ sigma^2;
      channelInput = (1 + exp(-2*y(j)*1/sigma^2))^(-1);% P(xj=1|yi)
      Lj(j) = channelInput;
      Prob = abs(channelInput/Nmax);
      initVNj(j,:) = (rand(1,LBernoulli)<Prob);       
    end
    
        % Initialiation    
    Z_ItoA(1,:) = initVNj(1,:);
    Z_ItoA(2,:) = initVNj(2,:);
    Z_ItoA(3,:) = initVNj(3,:);    
    Z_ItoA(4,:) = initVNj(3,:);
    Z_ItoA(5,:) = initVNj(4,:);
    Z_ItoA(6,:) = initVNj(4,:);
    Z_ItoA(7,:) = initVNj(4,:);
    Z_ItoA(8,:) = initVNj(5,:);
    Z_ItoA(9,:) = initVNj(5,:);
    Z_ItoA(10,:) = initVNj(6,:);
    Z_ItoA(11,:) = initVNj(6,:);
    Z_ItoA(12,:) = initVNj(7,:); 
    
    for BPiter = 1:BPiterMax
      %2. a) CN to VN updates
      Z_AtoI(1,:) = mod(Z_ItoA(3,:) + mod(Z_ItoA(5,:)+Z_ItoA(8,:),2), 2 ); %3,5,8
      Z_AtoI(2,:) = mod(Z_ItoA(2,:) + mod(Z_ItoA(5,:)+Z_ItoA(8,:),2), 2 ); %2,5,8
      Z_AtoI(3,:) = mod(Z_ItoA(2,:) + mod(Z_ItoA(3,:)+Z_ItoA(8,:),2), 2 ); %2,3,8
      Z_AtoI(4,:) = mod(Z_ItoA(2,:) + mod(Z_ItoA(3,:)+Z_ItoA(5,:),2), 2 ); %2,3,5
    
      Z_AtoI(5,:) = mod(Z_ItoA(6,:) + mod(Z_ItoA(9,:)+Z_ItoA(10,:),2), 2 ); %6,9,10
      Z_AtoI(6,:) = mod(Z_ItoA(1,:) + mod(Z_ItoA(9,:)+Z_ItoA(10,:),2), 2 ); %1,9,10
      Z_AtoI(7,:) = mod(Z_ItoA(1,:) + mod(Z_ItoA(6,:)+Z_ItoA(10,:),2), 2 ); %1,6,10
      Z_AtoI(8,:) = mod(Z_ItoA(1,:) + mod(Z_ItoA(6,:)+Z_ItoA(9,:),2), 2 ); %1,6,9
      Z_AtoI(9,:) = mod(Z_ItoA(7,:) + mod(Z_ItoA(11,:)+Z_ItoA(12,:),2), 2 ); %7,11,12
      Z_AtoI(10,:) = mod(Z_ItoA(4,:) + mod(Z_ItoA(11,:)+Z_ItoA(12,:),2), 2 ); %4,11,12
      Z_AtoI(11,:) = mod(Z_ItoA(4,:) + mod(Z_ItoA(7,:)+Z_ItoA(12,:),2), 2 ); %4,7,12
      Z_AtoI(12,:) = mod(Z_ItoA(4,:) + mod(Z_ItoA(7,:)+Z_ItoA(11,:),2), 2 ); %4,7,11      
        
      %2. b) VN to CN updates
      %     part 1 - Construct W
      %EQ_Z_i = EQopChan(initVNj(:,:));
      W_ItoA(1,:) = initVNj(1,:);
      W_ItoA(2,:) = initVNj(2,:);
      W_ItoA(3,:) = EQop(Z_AtoI(9,:), initVNj(3,:));  
      W_ItoA(4,:) = EQop(Z_AtoI(2,:), initVNj(3,:));
      W_ItoA(5,:) = EQop3(Z_AtoI(6,:),Z_AtoI(10,:), initVNj(4,:));
      W_ItoA(6,:) = EQop3(Z_AtoI(3,:),Z_AtoI(10,:), initVNj(4,:));
      W_ItoA(7,:) = EQop3(Z_AtoI(3,:),Z_AtoI(6,:), initVNj(4,:));
      W_ItoA(8,:) = EQop(Z_AtoI(7,:), initVNj(5,:));
      W_ItoA(9,:) = EQop(Z_AtoI(4,:), initVNj(5,:));
      W_ItoA(10,:) = EQop(Z_AtoI(11,:), initVNj(6,:));
      W_ItoA(11,:) = EQop(Z_AtoI(8,:), initVNj(6,:));
      W_ItoA(12,:) = initVNj(7,:);
    
      %   part 2 - Draw IID from W
      k = LBernoulli/2;
      IID_Samples = int32(rand(1,2*k)*k);
      Z_ItoA(1,:) = W_ItoA(1,[IID_Samples]+k);IID_Samples = int32(rand(1,2*k)*k);
      Z_ItoA(2,:) = W_ItoA(2,[IID_Samples]+k);IID_Samples = int32(rand(1,2*k)*k);
      Z_ItoA(3,:) = W_ItoA(3,[IID_Samples]+k);IID_Samples = int32(rand(1,2*k)*k);
      Z_ItoA(4,:) = W_ItoA(4,[IID_Samples]+k);IID_Samples = int32(rand(1,2*k)*k);
      Z_ItoA(5,:) = W_ItoA(5,[IID_Samples]+k);IID_Samples = int32(rand(1,2*k)*k);
      Z_ItoA(6,:) = W_ItoA(6,[IID_Samples]+k);IID_Samples = int32(rand(1,2*k)*k);
      Z_ItoA(7,:) = W_ItoA(7,[IID_Samples]+k);IID_Samples = int32(rand(1,2*k)*k);
      Z_ItoA(8,:) = W_ItoA(8,[IID_Samples]+k);IID_Samples = int32(rand(1,2*k)*k);
      Z_ItoA(9,:) = W_ItoA(9,[IID_Samples]+k);IID_Samples = int32(rand(1,2*k)*k);
      Z_ItoA(10,:) = W_ItoA(10,[IID_Samples]+k);IID_Samples = int32(rand(1,2*k)*k); 
      Z_ItoA(11,:) = W_ItoA(11,[IID_Samples]+k);IID_Samples = int32(rand(1,2*k)*k);
      Z_ItoA(12,:) = W_ItoA(12,[IID_Samples]+k); 
      
    end
    
    % c - compute marginal estimates
    U_i(1,:) = EQop(Z_AtoI(5,:), initVNj(1,:));
    U_i(2,:) = EQop(Z_AtoI(1,:), initVNj(2,:));
    U_i(3,:) = EQop3(Z_AtoI(2,:),Z_AtoI(9,:), initVNj(3,:));
    U_i(4,:) = EQop4(Z_AtoI(3,:),Z_AtoI(6,:),Z_AtoI(10,:), initVNj(4,:));
    U_i(5,:) = EQop3(Z_AtoI(4,:),Z_AtoI(7,:), initVNj(5,:));
    U_i(6,:) = EQop3(Z_AtoI(8,:),Z_AtoI(11,:),initVNj(6,:));
    U_i(7,:) = EQop(Z_AtoI(12,:),initVNj(7,:));
    
    n_i(1) = (1/k)*sum(U_i(1,(k+1:2*k)));
    n_i(2) = (1/k)*sum(U_i(2,(k+1:2*k)));
    n_i(3) = (1/k)*sum(U_i(3,(k+1:2*k)));
    n_i(4) = (1/k)*sum(U_i(4,(k+1:2*k)));
    n_i(5) = (1/k)*sum(U_i(5,(k+1:2*k)));
    n_i(6) = (1/k)*sum(U_i(6,(k+1:2*k)));
    n_i(7) = (1/k)*sum(U_i(7,(k+1:2*k)));   
    
    v_hat = n_i < 0.5;
    BitErrs = length(find(zeros(1,7) ~= v_hat(1:end)));
    
    num = num + L;
    totErr = totErr + BitErrs;
    disp(['iteration finished, totErr: ' num2str(totErr) ' num: ' num2str(num)]);
  end
  BER(i) = totErr/num;
end
%Uncoded Transmission
uncoded = qfunc(sqrt(2*EbN0_lin));
 
figure
semilogy(EbN0,uncoded)
hold on;
semilogy(EbN0, BER)
str = sprintf('Coded - Simulation. Iterations: %d',maxIter);

BER