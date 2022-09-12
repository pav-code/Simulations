%% Transmitter %%

function [SM,z] = CodedTransmitter(symbols, N, Nsym, Ncp, Ts, diversity)

% alphabet = 65:(65+26);
% for i = 1:length(alphabet)
%    ascii(i,:) = dec2bin(alphabet(i)); 
% end

SM = reshape(symbols,N,[]).'; % The symbols are devided to the sub-carriers

OFDM = sqrt(N/Ts).*ifft(SM,[],2);  % We apply inverse-fft on our symbols (OFDM)

Prefix = OFDM(:, end-Ncp+1:end); % Cyclic prefix -> Enables the properties for circular convolution
                              
z_long = [Prefix OFDM]; % We add the prefix to the signal

z = reshape(z_long.',Nsym*diversity*(N+Ncp),[]); % We reshape the signal into one long vector

end