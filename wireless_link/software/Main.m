function [eyed,const,d] = main(N)
close all;
clear all;
c = config;

% Excluded:
%  factor : [receiverCheckBarker.m]

% ====  General Variables ==== 
msgLength        = 1600;
freqOffset       = 700e3 + 100000*rand(1);
phaseOffset      = 2*pi*rand(1);%(3*pi/4 + 100*rand(1)); 0.2;
pilotFreq        = 2e6;
noise            = 15;
fc               = 4e6;
samp_f           = c.samplingFreq;
sym_rate         = c.symbolRate;
alpha_rollOff    = c.alpha; %(0.8)
rrtSpan          = c.span; %(8)
samplesPerSymbol = ceil(samp_f/sym_rate);
samPerSymRaw     = samp_f/sym_rate;
mapping          = 'QPSK';
barker_seq       = c.barker;

% ====  Rx Specific ==== 
N                = 1048576;

%% Create Bits
bits = randsrc(1,msgLength-60,[0 1]);
x = [0 0 1 1];
x = repmat(x,1,15);
bits = [x bits];                   %Add 60 bits of alternating 0011 quartet 


%% Transmitter
symbols = bits2symbols(bits,mapping);          %0s,1s to QPSK

symbolsI = real(symbols);                      %To 1s,-1s
symbolsQ = imag(symbols);

symbolsI = [c.barker c.barker symbolsI];       %Add Barker
symbolsQ = [c.barker c.barker symbolsQ];

symbolsI = awgn(symbolsI,noise);               %Add Noise
symbolsQ = awgn(symbolsQ,noise);

[rtrcPulse,~] = rtrcpuls(alpha_rollOff,1/sym_rate,samp_f,rrtSpan); %Grab a RRC Pulse

upsampledSymbolsCos = upsample(symbolsI,samplesPerSymbol); %Upsample to samples/symbols
I = conv(upsampledSymbolsCos, rtrcPulse);                  %Layout the RRCs

upsampledSymbolsSin = upsample(symbolsQ,samplesPerSymbol); %Quadrature part
Q = conv(upsampledSymbolsSin, rtrcPulse);

t_vec = 0: 1/samp_f: ((length(I)/samp_f)-(1/samp_f));      %Modulate to carrier freq
ptCos = I.*cos(2*pi*(fc)*t_vec+ phaseOffset);
ptSin = Q.*sin(2*pi*(fc)*t_vec+ phaseOffset);

pilotSin = 0.07*sin(2*pi*(fc+pilotFreq)*t_vec);            %Crude freq estimator
%Caution "2*10" mystery factor
passband = sqrt(2*10*sym_rate)*(ptCos+ptSin+pilotSin); %Signals. 


%% Channel
passband = passband.*exp(1j*2*pi*(freqOffset)*t_vec);

%% Receiver

%% 1) Crude Frequency Estimation
pilot = fft(passband,N);                                %Retreive pilot freq offset
[freqMax,freqInd] = max(abs(pilot(1:ceil(end/2))));         %Find the peak, and frequency
initFreqoffset = (freqInd-1)*(samp_f/N);                    %Recover freq. from FFT index
initFreqoffsetRel = initFreqoffset-(fc+pilotFreq);            %Get relative freq shift
t1 = 0:1/samp_f:(length(passband)/samp_f)-(1/samp_f);       %Repare for correction
passband_st1 = passband.*exp(1j*2*pi*(-initFreqoffsetRel)*t1); %Correction!


t = 0:samp_f/N:samp_f-1/(samp_f/N);         %Graphs...
figure
plot(t,abs(pilot)); hold on;
plot(initFreqoffset,abs(freqMax),'r*'); hold off;
xlim([0 c.samplingFreq]);
title('Pilot signal');

p2 = fft(passband_st1,N);                   %More graphs
figure
plot(t,abs(p2));
xlim([0 c.samplingFreq/2]);
title('Baseband pilot');
    
disp('Crude frequency offset')                              %Outputs to prompt
disp(initFreqoffset)
disp('Crude difference from Actual')
disp(initFreqoffsetRel - freqOffset)

%% Filter
stopShift = 2.5e5; passShift = 7.5e5;                       %Frequencies stop and pass shift limits in frequency

passNorm = (fc+stopShift)/(samp_f/2);                       %Normalized frequencies for filter
stopNorm = (fc+passShift)/(samp_f/2);                       %stop and pass bands
lpFilt = designfilt('lowpassfir','PassbandFrequency',passNorm, ... %Generating the filter
    'StopbandFrequency',stopNorm,'PassbandRipple',0.5, ...
    'StopbandAttenuation',40,'DesignMethod','kaiserwin');
passband_st1 = filter(lpFilt,passband_st1);                 %Filtering!

p3 = fft(passband_st1,N);                               %Graphs...
figure
plot(t,abs(p3));
xlim([0 c.samplingFreq/2]);
title('Baseband pilot filtered');

disp('Passband Edge')                                   %Outputs to prompt
disp(passNorm)
disp('Stopband Edge')
disp(stopNorm)

%% 2) Fine frequency correction
cosfc = cos(2*pi*(fc)*t1);                          %Back to pulse shapes for Bark autocorrelation
sinfc = sin(2*pi*(fc)*t1);
sI_baseband = passband_st1.*cosfc;
sQ_baseband = passband_st1.*sinfc;
sI_mfr = conv(sI_baseband,rtrcPulse);
sQ_mfr = conv(sQ_baseband,rtrcPulse);

[barker,~] = symb2rrc(c.barker,c.barker);           %construct RRC Barker
barker = [barker barker];                           %Double barker
sI_corr = conv(barker,fliplr(sI_mfr));              %Convolution for autocorrelation
sQ_corr = conv(barker,fliplr(sQ_mfr));

corr = sqrt(sI_corr.^2 + sQ_corr.^2);               %Autocorrelate
[~,index] = max(abs(corr));                         %Find synch point

figure                                              %Graph...
plot(abs(corr));hold on;
plot(index,max(abs(corr)),'r*'); hold off;

framestart = length(sI_mfr)-ceil((index-(29.99*samPerSymRaw))); %Frame starting found.
sI_frame = sI_mfr(framestart-ceil(22*samPerSymRaw):end);        %Gather all the symbol points
sQ_frame = sQ_mfr(framestart-ceil(22*samPerSymRaw):end);

sampl = (((sI_frame(1+(0:799)*samPerSymRaw)+ ...                %...and see their drift (all are folded onto angle pi)
    (1j*sQ_frame(1+(0:799)*samPerSymRaw)))).*exp(1j*pi/4)).^4;
angl = unwrap(angle(sampl));
x = 1:length(angl);                                 %Fitting a line of best fit through: angle vs sample number
p = polyfit(x,angl,1);
freqoffset = p(1)/4;                                %Divide by 4 since we raised to the power of 4, previously.

passband_st2 = passband_st1.*exp(1j*2*pi*(-freqoffset * samPerSymRaw * 100)*t1); %Make the fine correction!
disp('Fine Correction')
disp(-freqoffset * samPerSymRaw * 100)



%% 3) Phase correction
% Do the same thing as 2), to get the new index point and do the final
% correction, on passband_st2: passband stage 2.
cosfc = cos(2*pi*(fc)*t1);                          %Back to pulse shapes for Bark autocorrelation
sinfc = sin(2*pi*(fc)*t1);
sI_baseband = passband_st2.*cosfc;
sQ_baseband = passband_st2.*sinfc;
sI_mfr = conv(sI_baseband,rtrcPulse);
sQ_mfr = conv(sQ_baseband,rtrcPulse);
[barker,~] = symb2rrc(c.barker,c.barker);           %construct RRC Barker
barker = [barker barker];                           %Double barker
sI_corr = conv(barker,fliplr(sI_mfr));              %Convolution for autocorrelation
sQ_corr = conv(barker,fliplr(sQ_mfr));
corr = sqrt(sI_corr.^2 + sQ_corr.^2);               %Autocorrelate
[~,index] = max(abs(corr));                         %Find synch point
framestart = length(sI_mfr)-ceil((index-(29.99*samPerSymRaw))); %Frame starting found.
sI_frame = sI_mfr(framestart-ceil(22*samPerSymRaw):end);        %Gather all the symbol points
sQ_frame = sQ_mfr(framestart-ceil(22*samPerSymRaw):end);

phase = [angle(barker_seq + 1j*barker_seq)    ...           %Average the angles to find the phase
    angle(barker_seq + 1j*barker_seq)] - ...
    angle(sI_frame(1 + samPerSymRaw.*[0:25]) + ...
    1j*sQ_frame(1 + samPerSymRaw.*[0:25]));

for R = 1:length(phase)                             %Make negative phases positive (MATLAB annoyance)
    if sign(phase(R)) ~= 1 && phase(R) < -0.01
        phase(R) = 2*pi + phase(R);
    end  
end
phasedifference = sum(phase)/26;                    %CAUTION: We use double barker so this is 26 (CHANGE LATER)


disp('Phase Correction')
disp(phasedifference)


%% Final Stage: Recovery with corrected crude and fine frequency correction and Phase correction.
cosfc = cos(2*pi*(fc)*t1 + phasedifference);
sinfc = sin(2*pi*(fc)*t1 + phasedifference);

sI_baseband = passband_st2.*cosfc;
sQ_baseband = passband_st2.*sinfc;
sI_mfr = conv(sI_baseband,rtrcPulse);
sQ_mfr = conv(sQ_baseband,rtrcPulse);
sI_framed = sI_mfr(framestart-23*samPerSymRaw:end);%receiverFindFrame(sI_mfr);
sQ_framed = sQ_mfr(framestart-23*samPerSymRaw:end);%receiverFindFrame(sQ_mfr);

eyed.fsfd = samPerSymRaw;               %Eye Diagram
eyed.r = sI_framed(1:(c.messageLength/4 + 11)*samPerSymRaw) + 1i*sQ_framed(1:(c.messageLength/4 + 11)*samPerSymRaw);

sI_sample = downsample(sI_framed(1:end),samPerSymRaw); %TODO: sanity check below
sQ_sample = downsample(sQ_framed(1:end),samPerSymRaw);

sI_normalizationFactor = sum(abs(sI_sample(1:26)))/26;
sQ_normalizationFactor = sum(abs(sQ_sample(1:26)))/26;

sI_sample = sI_sample(27:c.messageLength/4+26);
sQ_sample = sQ_sample(27:c.messageLength/4+26);


%% Scale to correct size for constellation so ML receiver works
sI_symbols = (3*sI_sample)/(sqrt(10)*sI_normalizationFactor);
sQ_symbols = (3*sQ_sample)/(sqrt(10)*sQ_normalizationFactor);


s_bits = symbols2bits(sI_symbols+ 1j*sQ_symbols,mapping);

pack = s_bits;

psd = 0;
disp('Completed reception')

%% Calculate return values
const = sI_symbols+1i*sQ_symbols;

figure
plot(const,'r*')

eyediagram(eyed.r,eyed.fsfd)


isequal(bits,pack)
sum(abs(bits(1:40) - pack(1:40)))
bits(1,1:40)
end
