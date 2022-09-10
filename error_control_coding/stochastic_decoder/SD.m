%%% Decoder Stochastic %%%
clear all;
% Parameter
EbN0 = 0:1:3;
EbN0_lin = 10.^(EbN0/10);
R = 4/7;
L = 7; VNj = 7; CNi = 3;
maxNumErrs = 100;
maxNum = 1e5;

%SD parameters
Nmax = 1;
LBernoulli = 200;
BPiterMax = 15;

encTxF = [1 1 1 1 1 1 1]; %[all zero vector: BPSK] 

for i = 1:length(EbN0)
  totErr = 0;
  num = 0;
  while((totErr < maxNumErrs) && (num < maxNum)) 
    %Noise Addition

    sigma = 1/sqrt(2*EbN0_lin(i)*R);
    cgnoise = sigma.*randn(1,length(encTxF));
    %cgnoise = 0;
    y = encTxF + cgnoise;
      
    %Initializaton (1)
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
    E_EtoP(1,:) = initVNj(1,:);
    E_EtoP(2,:) = initVNj(2,:);
    E_EtoP(3,:) = initVNj(3,:);    
    E_EtoP(4,:) = initVNj(3,:);
    E_EtoP(5,:) = initVNj(4,:);
    E_EtoP(6,:) = initVNj(4,:);
    E_EtoP(7,:) = initVNj(4,:);
    E_EtoP(8,:) = initVNj(5,:);
    E_EtoP(9,:) = initVNj(5,:);
    E_EtoP(10,:) = initVNj(6,:);
    E_EtoP(11,:) = initVNj(6,:);
    E_EtoP(12,:) = initVNj(7,:); 
    
    for BPiter = 1:BPiterMax
    
    % Parity Check (PC) to equality nodes (EQ)
    E_PtoE(1,:) = PCNode(PCNode(E_EtoP(3,:),E_EtoP(5,:)),E_EtoP(8,:)); %3,5,8
    E_PtoE(2,:) = PCNode(PCNode(E_EtoP(2,:),E_EtoP(5,:)),E_EtoP(8,:)); %2,5,8
    E_PtoE(3,:) = PCNode(PCNode(E_EtoP(2,:),E_EtoP(3,:)),E_EtoP(8,:)); %2,3,8
    E_PtoE(4,:) = PCNode(PCNode(E_EtoP(2,:),E_EtoP(3,:)),E_EtoP(5,:)); %2,3,5
    E_PtoE(5,:) = PCNode(PCNode(E_EtoP(6,:),E_EtoP(9,:)),E_EtoP(10,:)); %6,9,10
    E_PtoE(6,:) = PCNode(PCNode(E_EtoP(1,:),E_EtoP(9,:)),E_EtoP(10,:)); %1,9,10
    E_PtoE(7,:) = PCNode(PCNode(E_EtoP(1,:),E_EtoP(6,:)),E_EtoP(10,:)); %1,6,10
    E_PtoE(8,:) = PCNode(PCNode(E_EtoP(1,:),E_EtoP(6,:)),E_EtoP(9,:)); %1,6,9
    E_PtoE(9,:) = PCNode(PCNode(E_EtoP(7,:),E_EtoP(11,:)),E_EtoP(12,:)); %7,11,12
    E_PtoE(10,:) = PCNode(PCNode(E_EtoP(4,:),E_EtoP(11,:)),E_EtoP(12,:)); %4,11,12
    E_PtoE(11,:) = PCNode(PCNode(E_EtoP(4,:),E_EtoP(7,:)),E_EtoP(12,:)); %4,7,12
    E_PtoE(12,:) = PCNode(PCNode(E_EtoP(4,:),E_EtoP(7,:)),E_EtoP(11,:)); %4,7,11

    % Equality Node (EQ) to Parity Check Nodes (PC) 
    
    E_EtoP(1,:) = initVNj(1,:);
    E_EtoP(2,:) = initVNj(2,:);
    E_EtoP(3,:) = EQNode(E_PtoE(9,:),initVNj(3,:));   
    E_EtoP(4,:) = EQNode(E_PtoE(2,:),initVNj(3,:)); 
    E_EtoP(5,:) = EQNode(EQNode(E_PtoE(6,:),E_PtoE(10,:)),initVNj(4,:));
    E_EtoP(6,:) = EQNode(EQNode(E_PtoE(3,:),E_PtoE(10,:)),initVNj(4,:));
    E_EtoP(7,:) = EQNode(EQNode(E_PtoE(3,:),E_PtoE(6,:)),initVNj(4,:));
    E_EtoP(8,:) = EQNode(E_PtoE(7,:),initVNj(5,:));
    E_EtoP(9,:) = EQNode(E_PtoE(4,:),initVNj(5,:));
    E_EtoP(10,:) = EQNode(E_PtoE(11,:),initVNj(6,:));
    E_EtoP(11,:) = EQNode(E_PtoE(8,:),initVNj(6,:));
    E_EtoP(12,:) = initVNj(7,:); 
    
    end
    % Total LLR summation
    for q = 1:12
       Pparity = length(find(E_PtoE(q,:)==1))/LBernoulli;
       Lij(q) = Pparity * Nmax;
    end
    %Lij = Lij/12
    LjTotal(1) = 0.5*(Lj(1) + Lij(5));
    LjTotal(2) = 0.5*(Lj(2) + Lij(1));
    LjTotal(3) = 0.5*(Lj(3) + Lij(2) * Lij(9));
    LjTotal(4) = 0.5*(Lj(4) + Lij(3) * Lij(6) * Lij(10)); 
    LjTotal(5) = 0.5*(Lj(5) + Lij(7) * Lij(7));
    LjTotal(6) = 0.5*(Lj(6) + Lij(8) * Lij(11));
    LjTotal(7) = 0.5*(Lj(7) + Lij(12));
 
    v_hat = floor(-0.5*sign(LjTotal)+ 0.5); 
    v_hat = LjTotal < 0.5;
    
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