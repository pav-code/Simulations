function plotDice()
 
close all; clear all;
level       = [ 0         3      5       7       9 ];
diceGM      = [0.94     0.91    0.90    0.85    0.83 ];
diceWM      = [0.88     0.85    0.80    0.75    0.90 ]; 
diceCSF     = [0.99     0.97    0.94    0.90    0.86 ];
 
figure; LW=2; set(gca,'FontSize', 14);
plot(level, diceCSF,'ro','LineWidth', 2); grid on; hold on;
plot(level, diceGM,'ks','LineWidth', 2);  
plot(level, diceWM,'bd','LineWidth', 2); 
 

xlabel('Noise level [%]');
ylabel('Dice index');
axis([-1 10 0.6 1.0]);
legend('CSF', 'Gray Matter', 'White Matter');

plot(level, diceCSF,'r:');  
plot(level, diceGM,'k:'); 
plot(level, diceWM,'b:');  
title('Dice index for T1 and T2, slice #91, INU=20%'); 

%     ax = axis; xlabel('Dice index');
%     axis([0 8 ax(3) ax(4)]);
    set(gca,'XTick',[ 0 3 5 7 9]);
%     set(gca,'XTickLabel',[ {' '} {'CSF'} {'GM'} {'WM'} {'FAT'} {'MUSCLE'} {'SKIN'} {'SKULL'} {' '}]);
 
return
