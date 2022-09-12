M = 16; % number of constellation points
k = log2(M); % number of bits in each constellation
% def5ining the real and imaginary PAM constellation
% for 16-QAM
alphaRe = [-(2*sqrt(M)/2-1):2:-1 1:2:2*sqrt(M)/2-1];
alphaIm = [-(2*sqrt(M)/2-1):2:-1 1:2:2*sqrt(M)/2-1];
% input - decimal equivalent of all combinations with b0b1b2b3
ip = [0:15];
ipBin = dec2bin(ip.'); % decimal to binary
% taking b0b1 for real
ipDecRe = bin2dec(ipBin(:,[1:k/2]));
ipGrayDecRe = bitxor(ipDecRe,floor(ipDecRe/2));
% taking b2b3 for imaginary
ipDecIm = bin2dec(ipBin(:,[k/2+1:k]));
ipGrayDecIm = bitxor(ipDecIm,floor(ipDecIm/2));
% mapping the Gray coded symbols into constellation
modRe = alphaRe(ipGrayDecRe+1);
modIm = alphaIm(ipGrayDecIm+1);
% complex constellation
mod = modRe + j*modIm;



x = (0:63)';
y = bin2gray(x,'qam',64);

z = gray2bin(y,'qam',64);
% 
% hMod = comm.RectangularQAMModulator;
% symbols = constellation(hMod);
hMod = comm.RectangularQAMModulator('ModulationOrder',64, ...
    'SymbolMapping','binary');
symbols = constellation(hMod);
scatterplot(symbols,1,0,'b*');

