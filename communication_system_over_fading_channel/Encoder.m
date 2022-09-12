function y = Encoder(randBits, diversity, E, m)

QPSK = [-1-1i; -1+1i; 1-1i; 1+1i].*sqrt(E/2); % QPSK constellation

bits  = randBits; % Bits Input - diversified bits (frame spacing, diveristy) 
    
GroupBits = buffer(bits,m)'; % We group the bits 2 by 2

codeword = bi2de(GroupBits,'left-msb')+1; % Assign each "group" to a decimal

symbols = QPSK(codeword); % Each number is assigned to a point of the constellation

%newBits = zeros(length(symbols)*diversity,1);
newBits = zeros(length(symbols),1);
if diversity == 1
    newBits = symbols;
else
    for k = 1:length(symbols)
        for kk = 1:diversity
            newBits(diversity*(k-1)+kk,1) = symbols(k);
        end
    end
end

y = newBits;
end