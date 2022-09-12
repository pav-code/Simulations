%% Transmitter %%

function [SM,z,bits] = Transmitter(m, N, Nsym, Nbits, Ncp, Ts, E)

% alphabet = 65:(65+26);
% for i = 1:length(alphabet)
%    ascii(i,:) = dec2bin(alphabet(i)); 
% end

QPSK = [-1-1i; -1+1i; 1-1i; 1+1i].*sqrt(E/2); % QPSK constellation

bits  = randi([0 1],Nbits,1); % Bits generation
    
GroupBits = buffer(bits,m)'; % We group the bits 2 by 2

codeword = bi2de(GroupBits,'left-msb')+1; % Assign each "group" to a decimal

symbols = QPSK(codeword); % Each number is assigned to a point of the constellation

SM = reshape(symbols,N,[]).'; % The symbols are devided to the sub-carriers

OFDM = sqrt(N/Ts).*ifft(SM,[],2);  % We apply inverse-fft on our symbols (OFDM)

Prefix = OFDM(:, end-Ncp+1:end); % Cyclic prefix -> Enables the properties for circular convolution
                              
z_long = [Prefix OFDM]; % We add the prefix to the signal

z = reshape(z_long.',Nsym*(N+Ncp),[]); % We reshape the signal into one long vector

end