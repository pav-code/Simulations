function [c] = EQNode(a,b)
%Equality Node Equation
J = a&b;
K = not(a) & not(b);
%(1-Pa)*Pb + (1-Pb)*Pa + (1-Pa)*(1-Pb)

JKout = zeros(1,length(J));
for i = 1:length(J)
  if     J(i) == 0 && K(i) == 0
    if i == 1
      JKout(1) = 0;  
    else
      JKout(i) = JKout(i-1);
    end
  elseif J(i) == 0 && K(i) == 1
    JKout(i) = 0;  
  elseif J(i) == 1 && K(i) == 0
    JKout(i) = 1;  
  elseif J(i) == 1 && K(i) == 1
    if i == 1
      JKout(1) = 1;  
    else
      JKout(i) = not(JKout(i-1));
    end      
  end
end
c = JKout;
end