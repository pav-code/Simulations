function y = DeintleaverTime(N, bits, diversity, frameSpacing, Nsym)

index = 1; fIndex = 1; outInx = 1; base = 0;
fSize = frameSpacing*diversity*N; %(fSpace + 1)*N
fSize = (frameSpacing + 1)*N;
y = zeros(length(bits),1);

for i = 1:diversity:length(bits)
  for j = 1:diversity
    
    offset = index + N*frameSpacing*(j-1) + (fIndex-1)*N + base;
    %y(offset) = bits(i);
    y(outInx) = bits(offset);
    outInx = outInx + 1;
      
  end
  if N == index
    index = 0;
    fIndex = fIndex + 1;
    if fIndex > frameSpacing
      fIndex = 1;
      base = base + diversity*frameSpacing*N;
    end
  end
  index = index + 1;
end

end