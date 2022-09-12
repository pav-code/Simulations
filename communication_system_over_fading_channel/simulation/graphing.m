d1 = load('BER_d1_graph.mat');
d1 = d1.BER;
d2 = load('BER_d2_graph.mat');
d2 = d2.BER;
d3 = load('BER_d3_graph.mat');
d3 = d3.BER;
EbN0 = [0:1:20];
EbN0lin = 10.^(EbN0/10);

xi=0:1:20;
y1=interp1(EbN0,d1,xi,'spline');

xi=0:1:20;
y2=interp1(EbN0,d2,xi,'spline');

xi=0:1:20;
y3=interp1(EbN0,d3,xi,'spline');

uncoded = qfunc(sqrt(2*EbN0lin)); % AWGN QPSK, BER curve (comparison)

figure
semilogy(EbN0,uncoded)
hold on;
semilogy(EbN0, y1, 'bo')
semilogy(EbN0, y2, 'ro')
semilogy(EbN0, y3, 'ko')
ylim([10^-6 1])
title('Uncoded Rayleigh fading channel')
ylabel('BER')
xlabel('Eb/N0 [dB]')
legend('AWGN', 'Rayleigh with 101 [dB] P_L', 'a', 'b')