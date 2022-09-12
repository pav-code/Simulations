%% Receiver %%
function y = SoftReceiver(r, h, N, Nsym, Ncp, Ts)

% Two conditions for equations (30):
%   1) (N + Nsym)*fdTs << 1
%   2) Ncp >= ceil( delay L, Tau )

r = r(1:(N+Ncp)*Nsym); % We remove the length added by the channel

r_matrix = reshape(r,N+Ncp,[]).'; % The symbols are devided to the sub-carriers

r_cut = r_matrix(:,Ncp+1:end); % We remove the cycle prefix

data = sqrt(Ts/N).*fft(r_cut,[],2); % We apply fft on our symbols (OFDM)

% Generatiion of the C^m for channel estimation
H = [h(1,1) 0 0 0 h(1,2)]; % adding zeros due to the spacing of the delays.
                           % i.e.: tau1=0, tau2=4. Therefore the channel
                           % impulse response is 0 in between.
C = diag( fft( H, N) );

y = zeros(Nsym,N);

y = transpose(C'*data.'); % We multiply the received vector with C hermitian

%y = data; % by-passing channel test

yy = reshape(y.',1,[]);% We reshape the received signal into one long vector

y = yy;
end