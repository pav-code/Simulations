clear all;
% P(A,B) = P(A)*P(B) [independence]
% P(A,B) = Za AND Zb

for i = 1:100
L=500;
Pa = 0.25;
Pb = 0.778;

%Generate Bernoulli sequences
Za = (rand(1,L)<Pa);
Zb = (rand(1,L)<Pb);
Zc = Za&Zb;

Pc = length(find(Zc==1))/L;

estimate = Pc;
actual = Pa * Pb;
avEstimate(i) = estimate
percentError(i) = abs(estimate-actual) * 100;
%i
end
averageError = sum(percentError)/100;
aveEst = sum(avEstimate)/100;

%Parity Check Equation
Zparity = bitxor(Za,Zb);
Pparity = length(find(Zparity==1))/L;
actualPar = Pa*(1-Pb) + Pb*(1-Pa);
pErrParity = abs(Pparity-actualPar)

%Equality Node Equation
J = Za&Zb;
K = not(Za) & not(Zb);
%(1-Pa)*Pb + (1-Pb)*Pa + (1-Pa)*(1-Pb)

JKout = zeros(1,length(J));
for i = 1:length(J)
  if     J(i) == 0 && K(i) == 0
    if i == 1
      JKout(1) = 0;  
    else
      JKout(i) = JKout(i-1);
    end
  elseif J(i) == 0 && K(i) == 1
    JKout(i) = 0;  
  elseif J(i) == 1 && K(i) == 0
    JKout(i) = 1;  
  elseif J(i) == 1 && K(i) == 1
    if i == 1
      JKout(1) = 1;  
    else
      JKout(i) = not(JKout(i-1));
    end      
  end
end

actualEq = (Pa*Pb) / (1 - Pa - Pb + 2*Pa*Pb);
Peq = length(find(JKout==1))/L;
pErrEquality = abs(Peq-actualEq)


%%% Decoder Stochastic %%%



%%% Decoder SPA %%%
maxIter = 20;
maxNumErrs = 300;
maxNum = 1e6;
L = 7;
R = 4/7;
encTxF = [1 1 1 1 1 1 1]; %[all zero vector: BPSK]  
CNi = 3;
VNj = 7;

LijIndex = [1 3 5 7; 2 3 6 7; 4 5 6 7];
LjiIndex = [1 0 0; 0 2 0; 1 2 0; 0 0 3; 1 0 3; 0 2 3; 1 2 3];
%Noise Bounds
EbN0 = 0:1:7;
%EbN0 = 1:1:4;
EbN0_lin = 10.^(EbN0/10);

%Simulation Loop
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
    initVNj = zeros(1,VNj);
    %AWGN
    for j = 1:1:VNj
      initVNj(j) = 2*y(j) ./ sigma^2;
    end
    Lji(1,1) = initVNj(1);
    Lji(2,2) = initVNj(2);
    Lji(3,1) = initVNj(3);
    Lji(3,2) = initVNj(3);
    Lji(4,3) = initVNj(4);
    Lji(5,1) = initVNj(5); 
    Lji(5,3) = initVNj(5); 
    Lji(6,2) = initVNj(6); 
    Lji(6,3) = initVNj(6); 
    Lji(7,1) = initVNj(7); 
    Lji(7,2) = initVNj(7); 
    Lji(7,3) = initVNj(7); 
    
    numIter = 0;
     while 1 %Iteration loop
      %Lij (2) - magnitude reprisents the certainty of the bit according to
      %the SPC constraint on CNi (excluding VNj)
      mul = 1;
      for CNindex = 1:CNi
        for VNc = 1:size(LijIndex,2)
          VNindex = LijIndex(CNindex, VNc);
          for jPrimeCount = 1:4
            jPrime = LijIndex(CNindex,jPrimeCount);
            neighbourCheck = find(LijIndex(CNindex,:) == jPrime);
            if jPrime == VNindex
              continue;
            elseif isempty(neighbourCheck)
                break
            end
            mul = mul* tanh(0.5*Lji(jPrime, CNindex));
          end
          if isempty(neighbourCheck) ~= 1
            Lij(CNindex,VNindex) = 2*atanh(mul);
            mul = 1;
          end
        end
      end
      
      %Lji (3) - magnitude represents the certainty of the bit. Based on
      %REP from VNj, excluding CNi.
      sum1 = 0;
      for VNindex = 1:VNj
        for CNc = 1:size(LjiIndex,2)
          CNindex = LjiIndex(VNindex, CNc);
          for iPrimeCount = 1:3
             iPrime = LjiIndex(VNindex, iPrimeCount);
             if iPrime == CNindex || iPrime == 0
               continue; 
             end
             sum1 = sum1 + Lij(iPrime, VNindex);
          end
          if CNindex ~= 0
            Lji(VNindex,CNc) = initVNj(VNindex) + sum1;
          end
          sum1 = 0;
        end
      end
      
      %Lj_total (4)
      for VNindex = 1:VNj
        for iPrime = 1:3
          CNindex = LjiIndex(VNindex,iPrime);
          if CNindex ~= 0
            sum1 = sum1 + Lij(CNindex, VNindex); 
          end          
        end
        LjTotal(VNindex) = initVNj(VNindex) + sum1;
        sum1 = 0;
      end
      
      %v_hat (5)
      v_hat = floor(-0.5*sign(LjTotal)+ 0.5);
      %v_hat = LjTotal > 0;
      
%       %stopping criterion (5)
%       stopCri = mod(v_hat*H_original',2);
%       if ZeroVec(1) == stopCri(1) ...
%       && ZeroVec(2) == stopCri(2) ...
%       && ZeroVec(3) == stopCri(3) ...
%       || numIter > maxIter
%          break;
%       end
      if numIter > maxIter
         break;
      end
      numIter = numIter + 1;
    end
    %v_hat = floor(-0.5*sign(y)+ 0.5); %uncoded testing
    BitErrs = length(find(zeros(1,7) ~= v_hat(1:end)));
    
    num = num + L;
    totErr = totErr + BitErrs;
    disp(['iteration finished, totErr: ' num2str(totErr) ' num: ' num2str(num)]);
  end
  BER(i)    = totErr/num;
end
%Uncoded Transmission
uncoded = qfunc(sqrt(2*EbN0_lin));
 
figure
semilogy(EbN0,uncoded)
hold on;
semilogy(EbN0, BER)
str = sprintf('Coded - Simulation. Iterations: %d',maxIter);


%%%
% % 4 7 Hamming - ML %
% G = [1 1 0 1;
%     1 0 1 1;
%     1 0 0 0;
%     0 1 1 1;
%     0 1 0 0;
%     0 0 1 0;
%     0 0 0 1];
% x = [0 0 0 0;
%      0 0 0 1;
%      0 0 1 0;
%      0 0 1 1;
%      0 1 0 0;
%      0 1 0 1;
%      0 1 1 0;
%      0 1 1 1;
%      1 0 0 0;
%      1 0 0 1
%      1 0 1 0;
%      1 0 1 1;
%      1 1 0 0;
%      1 1 0 1;
%      1 1 1 0;
%      1 1 1 1];
%  C = mod(x*(G'),2);