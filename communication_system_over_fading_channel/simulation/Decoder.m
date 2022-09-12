function y = Decoder(bits, diversity)

% % outBits = zeros(length(bits)/diversity,1);
% % if (diversity == 1)
% %    outBits = bits;
% % else
% %     for k = 1:ceil(length(bits)/diversity)
% %         for kk = 1:diversity
% %             tmp(kk)=bits(diversity*k-kk+1,1);
% %         end
% %         outBits(k) = sum(tmp);
% %     end
% % end
decTot = 0;
for i = 1:diversity
  decTot = decTot + bits(i:diversity:end);
end

yy = decTot;

y_real = real(yy);
y_imag = imag(yy);

bits = y_real > 0;
bits2 = y_imag > 0;

bits3 = 2*ones(length(bits)*2,1);
bits3(1:2:2*length(bits)) = bits;
bits3(2:2:2*length(bits)) = bits2;

y = bits3;
end