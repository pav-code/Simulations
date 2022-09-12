%clear all;
%close all;
%clc;
%% Main Program %% 

%%%%% PARAMETERS %%%%%
diversity = 1;          % rate (n,k) code
frameSpacing = 2;       % frame spacing.
% Number of bits / (2*frameSpacing) must be an integer!!
% frameSpacing*N = Nbits 
fc   = 2e9;             % Carrier frequency   
fs   = 1e6;             % Sample frequency
Ts   = 1./fs;           % Symbol period (size of a QPSK symbol)
fd   = 15/(3e8/fc);     % Doppler frequency (with v = 15 m/s)
fdTs = fd*Ts;
Pt   = 0.1;             % Transmit power
Pr   = 7.94e-12;        % Receivered power after channel
m    = 2;               % Bits per symbol (QPSK)
N    = 128;             % Number of subcarriers (arbitrary chosen)
                        % Must be power of 2 (due to FFT speed up)
Nsym = 10000;           % Number of OFDM symbols to transmit
Nbits = 2 * N * Nsym;   % Total Number of bits
tau  = [0; 4];          % Channel delays
Ncp = ceil(max(tau))+6; % Length of the c.prefix ( > length of 
                        % channel delay spread (4 microseconds))
E = Ts*Pt;              % Transmitted energy per second joule
Tb = 5.405e-7;          % Tb effective. Tb = 1/Rb. Rb = Rb,eff * (N + Ncp) / N 
sigma = sqrt((2.07e-20*2)/Ts); % noise spectral density receiver, from PM (part 1)

EbN0 = [0:1:20];       % [Eb/N0] range of 0 - 25 dB. (SNR)
EbN0lin = 10.^(EbN0/10); % converning SNR dB to linear SNR

bits  = randi([0 1],Nbits,1); % Bits generation
%bits = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
%21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40]; %for debugging!
%%%%% ENCODER %%%%%
encSym = Encoder(bits, diversity, E, m); % Encoded bits

%%%%% INTERLEAVER %%%%%
interSym = InterleaverTime(N, encSym, diversity, frameSpacing, Nsym);
% BER Loop: run through different sigmas.
for i = 1:length(EbN0)
    
%%%%% NOISE %%%%%
N0 = Pr*Tb/EbN0lin(i);    %Part 2: BER vs EbN0 curve. 1) SNR = Eb/N0
                          % 2) N0 = Eb/SNR 3) N0 = Pr*Tb/SNR (SNR in linear
                          % scale)
sigma = sqrt(N0 / (2*Ts)); % w = white noise. PM p10, below (25). 
                           % variance(w) = N0/Ts. Sigma = sqrt(variance).
    
%%%%% TRANSMISSION %%%%%
[CM,z] = CodedTransmitter(interSym ,N, Nsym, Ncp, Ts, diversity);

Rx = zeros( length(interSym),1);
% Loop for sending the data through the channel OFDM symb per OFDM symb
for Tx = 1:Nsym*diversity
  currSymbol = z( (Tx-1)*(N+Ncp)+1: Tx*(N+Ncp) ); % Extract OFDM Symbol
 
  %%%%%% CHANNEL %%%%%%%
  %{
  Inputs:
  s = samples of the transmitted signal, M x 1
  tau = delays of taps in samples, L x 1
  fdTs = normalized Doppler frequency, Ts = sample interval
  P = power delay profile, optional, default is equal power taps

  % Outputs:
  r = noise-less output of channel, (M + max(tau)) x 1
  h = channel tap processes, (M + max(tau)) x L
  %}
  P = [0.5 0.5]; % Equal power distribution
  [r,h] = Fading_Channel(currSymbol, tau, fdTs); % Pass through channel

  r = r ./ sqrt(10^10.1); % add path loss - 101 [dB].
  
  %%%%% NOISE %%%%%
  w = sigma .* (randn(length(r),1) + 1i*randn(length(r),1)); % Generate AWGN
  r = r + w; % Add the noise to the post-channel signal.

  %%%%% RECEIVER %%%%%
  y = SoftReceiver(r, h, N, 1, Ncp, Ts);
  
  Rx( (Tx-1)*(N)+1: Tx*(N) ) = y; % We rebuild the received data
end
%%%%% DINTERLEAVER %%%%%%
deintBits = DeintleaverTime(N, Rx, diversity, frameSpacing, Nsym);

%%%%% DECODER %%%%%
Rx = Decoder(deintBits, diversity);

%%%%% SIMULATION %%%%%
percentErrors = sum(abs(Rx - bits))/length(bits);

disp(['SNR = ' num2str(EbN0(i)) ' [dB]. BER = ' num2str(percentErrors) ]);
BER(i) = percentErrors;
end

uncoded = qfunc(sqrt(2*EbN0lin)); % AWGN QPSK, BER curve (comparison)

figure
semilogy(EbN0,uncoded)
hold on;
grid on
title('Uncoded Rayleigh fading channel')
ylabel('BER')
xlabel('Eb/N0 [dB]')
legend('AWGN', 'Rayleigh with 101 [dB] P_L')

