function [DI,JI] = diceIndex(J)
% Dice index from J 2-d matrix
 
J = single(J);
[M N K] = size(J);

if (K>1)
    errmsg('diceIndex: the input matrix is more than 2D');
elseif (M~=N)
    errmsg('diceIndex: the input matrix is not quadratic');
elseif (M>0)
    for i=1:M
        D(i) = 2*J(i,i)/(  sum(J(i,:)) + sum(J(:,i)) );
        JI(i) = J(i,i)/(sum(J(i,:)) + sum(J(:,i)) - J(i,i));
    end   
end
 
if 1 % just for brain tissues
    figure; LW=2; set(gca,'FontSize', 14);
    plot(D(2:end) ,'ro:', 'LineWidth', LW); hold on; grid on; 
    %bar([1 2 3], D(2:end));
    ax = axis; xlabel('Tissue type'); ylabel('Dice index');
    axis([0 4 ax(3) ax(4)]);
    set(gca,'XTick',[ 0:4]);
    set(gca,'XTickLabel',[ {''} {'CSF'} {'GM'} {'WM'} {''}]);
    
    
end
DI = D;
return
