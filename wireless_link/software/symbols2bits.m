%% bits = symbols2bits(symbols,mapping)
%  Takes symbols and maps back into bits
%  given a mapping
function bits = symbols2bits(symbols,mapping)
    
    c = config();
    
    if strcmpi(mapping,'BPSK')
       bits = zeros(1,length(symbols));
       const = c.BPSK;
    elseif strcmpi(mapping,'QPSK')
       bits = zeros(2,length(symbols));
       const = c.QPSK;
    elseif strcmpi(mapping,'AMPM')
       bits = zeros(3,length(symbols));
       const = c.AMPM;
    elseif strcmpi(mapping,'16-QAM')
       bits = zeros(4,length(symbols));
       const = c.QAM;
    else
        disp('No/or unavailable mapping specified, please retry');
        bits = 0;
    end


   if strcmpi(mapping,'BPSK')
       indexes = zeros(length(symbols),1);
       for i = 1:size(bits,2)
           [~,symbol] = min(abs(symbols(i)-const));
           indexes(i) = symbol;
       end
          s = de2bi((indexes-1),1,2,'left-msb');
          bits = s';    
   elseif strcmpi(mapping,'QPSK')
%        indexes = zeros(length(symbols),1);
%        for i = 1:size(bits,2)
%            [~,symbol] = min(abs(symbols(i)-const));
%            indexes(i) = symbol;
%        end
       [~,indexes] =  min(abs(bsxfun(@minus,symbols.',const)),[],2);
       s = de2bi((indexes-1),2,2,'left-msb');
       %bits(:,i) = s'; 
       bits = reshape(s',1,length(symbols)*2);
   elseif strcmpi(mapping,'AMPM')
       indexes = zeros(length(symbols),1);
       for i = 1:size(bits,2)
           [~,symbol] = min(abs(symbols(i)-const));
           indexes(i) = symbol;
       end
       s = de2bi((indexes-1),3,2,'left-msb');
       bits = reshape(s',1,length(symbols)*3);
   elseif strcmpi(mapping,'16-QAM')
       [~,indexes] = min(abs(bsxfun(@minus,symbols.',const)),[],2);
       s = de2bi((indexes-1),4,2,'left-msb');
       bits = reshape(s',1,length(symbols)*4);
   end
   
    
%     if strcmpi(mapping,'BPSK')
%            return;
%     elseif strcmpi(mapping,'QPSK')
%            bits = reshape(bits,1,length(symbols)*2);
%     elseif strcmpi(mapping,'AMPM')
%            bits = reshape(bits,1,length(symbols)*3);
%     end


end