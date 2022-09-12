function [I,Q] = symb2rrc(symbolsI,symbolsQ)
    c = config();
    samplesPerSymbol = ceil(c.samplingFreq/c.symbolRate);
     
    % generate pulse shape
    [rtrcPulse,~] = rtrcpuls(c.alpha,1/c.symbolRate,c.samplingFreq,c.span);
    % rtrcPulse = rtrcPulse*sqrt(1/c.symbolRate);
    
    % I(t) = C_0 SUM_n a_n Sqrt(T) v(t-nT)
    upsampledSymbolsCos = upsample(symbolsI,samplesPerSymbol);
    I = conv(upsampledSymbolsCos, rtrcPulse);
    
    % Q(t) = C_0 SUM_n a_n Sqrt(T) v(t-nT)
    upsampledSymbolsSin = upsample(symbolsQ,samplesPerSymbol);
    Q = conv(upsampledSymbolsSin, rtrcPulse);
end