
%% Config file
%  This file contains all Static Global variables used

function c = config
% Sampling frequency of USRP
    c.samplingFreq = 100e6;

% Parameters for root raised cosine
    c.bandwidth = 500000;
    c.span = 8;
    c.symbolRate = 250000;
    c.alpha = 0.8;
    c.shape = 'sqrt';

% Misc. Parameters
    c.messageLength = 1600;
% Barker Sequence for signal detection
c.barker = [ +3 +3 +3 +3 +3 -3 -3 +3 +3 -3 +3 -3 +3];
c.barker = c.barker/mean(abs(c.barker));

% QPSK Mapping
% 01         11
%
%
% 00         10
% Symbol Mapping
% 00 -> 0 ->  -1-i1
% 01 -> 1 ->  -1+i1
% 10 -> 2 ->  1-i1
% 11 -> 3 ->  1+i1
c.QPSK = [-1-1i -1+1i 1-1i 1+1i];
%c.QPSK = c.QPSK/abs(1+1i);   

%16-QAM Mapping
% 0000   0100    1100    1000
%
% 0001   0101    1101    1001
%
% 0011   0111    1111    1011
%
% 0010   0110    1110    1010
% 
% Symbol Mapping:
% 0000 -> 0  -> -3+3j
% 0001 -> 1  -> -3+1j
% 0010 -> 2  -> -3-3j
% 0011 -> 3  -> -3-1j
% 0100 -> 4  -> -1+3j
% 0101 -> 5  -> -1+1j
% 0110 -> 6  -> -1-3j
% 0111 -> 7  -> -1-1j
% 1000 -> 8  ->  3+3j
% 1001 -> 9  ->  3+1j
% 1010 -> 10 ->  3-3j
% 1011 -> 11 ->  3-1j
% 1100 -> 12 ->  1+3j
% 1101 -> 13 ->  1+1j
% 1110 -> 14 ->  1-3j
% 1111 -> 15 ->  1-1j
b = 1; 
c.QAM = [(-3*b)+(3j*b) (-3*b)+(1j*b) (-3*b)-(3j*b) (-3*b)-(1j*b) ...
    (-1*b)+(3j*b) (-1*b)+(1j*b) (-1*b)-(3j*b) (-1*b)-(1j*b) ...
    (3*b)+(3j*b) (3*b)+(1j*b) (3*b)-(3j*b) (3*b)-(1j*b) ...
    (1*b)+(3j*b) (1*b)+(1j*b) (1*b)-(3j*b) (1*b)-(1j*b)];
c.QAM = c.QAM/mean(abs(c.QAM));

end