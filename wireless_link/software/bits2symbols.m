%% [Symbols] = bits2symbols(bits,'mapping');
%  This function takes a string of bits and maps
%  onto a pre-defined symbol map. Valid arguments
%  for 'mapping' are 'BPSK','QPSK','AMPM','16-QAM'.

function symbols = bits2symbols(bits,mapping)
c = config();
if strcmpi(mapping,'BPSK')
   %tmp = bi2de(bits','left-msb')+1;
   tmp = bits+1;
   symbols = c.BPSK(tmp);
elseif strcmpi(mapping,'QPSK')
   buf = buffer(bits,2);
   tmp = bi2de(buf','left-msb')+1;
   symbols = c.QPSK(tmp);     
elseif strcmpi(mapping,'AMPM')
   buf = buffer(bits,3);
   tmp = bi2de(buf','left-msb')+1;
   symbols = c.AMPM(tmp);
elseif strcmpi(mapping,'16-QAM')
    buf = buffer(bits,4);
    tmp = bi2de(buf','left-msb')+1;
    symbols = c.QAM(tmp);
else
    disp('No/or unavailable mapping specified, please retry');
    symbols = 0;
end



end