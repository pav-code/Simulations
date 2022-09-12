clear all;

M = 3;
SNR = 0;
gamma = sqrt(SNR/(1+SNR));
totM = 0;
thresold = 5*10^-3;
index = 0;

while 1
  if SNR > 20
      break
  end
  
  index = index + 1;
  EbN0(index) = SNR;
  SNRlin = 10^(SNR/10);
  gamma = sqrt(SNRlin/(1+SNRlin));
  A = ((1-gamma)/2)^M;
  for i = 1:M
    m = i-1;
    totM = totM + nchoosek(M-1+m,m)*((1+gamma)/2)^m;
  end
  Pb(index) = A*totM;
%   if abs(Pb - thresold) < 0.01
%     break 
%   end
  totM = 0;
  SNR = SNR + 0.001;

end

semilogy(EbN0,Pb,'color',[.25 .25 .25])
hold on
ylim([10^-6 1])
legend('Simulation - Diversity 1', ...
    'Simulation - Diversity 2',...
    'Simulation - Diversity 3',...
    'Theory - Diversity 1',...
    'Theory - Diversity 2',...
    'Theory - Diversity 3')
xlabel('Eb/N0');
ylabel('BER');
title('Diversity - Simulation and Theory');
