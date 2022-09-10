function Y = EQop3(X1, X2, X3)
   Y0 = 0.5<rand(1);
   
   for i = 1:length(X1)
     if X1(i) == 1 && X2(i) == 1 && X3(i) == 1
       Y(i) = 1;
     elseif X1(i) == 0 && X2(i) == 0 && X3(i) == 0
       Y(i) = 0;         
     else
       if i == 1
         Y(i) = Y0;
       else
         Y(i) = Y(i-1); 
       end
     end
   end

   
%      for i = 1:length(varargin{1})
%      for j = 1:nVarargs 
%        if varargin{1}()
%      end
%    end 
   
%       fprintf('Total number of inputs = %d\n',nargin);
%    
%    nVarargs = length(varargin);
%    fprintf('Inputs in varargin(%d):\n',nVarargs)
%    for k = 1:nVarargs
%       fprintf('   %d\n', varargin{k})
%    end
   