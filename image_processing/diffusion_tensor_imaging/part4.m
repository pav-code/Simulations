clear all; close all; clc;

load('data.mat')

figure
hold on;
plot(data(:,1),data(:,2))
plot(data(:,1),data(:,3))
plot(data(:,1),data(:,4))
legend('S Normal','S Tumor','S Noise')
title('Given data VS b-factor')

% TODO: add in the rest of the data and plot accordingly.

% 1) outliers can be removed

% 2) least-squares fit, using S(b)=S0exp(-bD) for b below 1000
% use a linear function of type: y = mx + b
% take the logarithm of the data as stated in the hint, F = linFun(x,data).
x0 = [0 0];
xdata = data(1:5,1);
ydata1 = log(data(1:5,2));
y_hat1 = lsqcurvefit(@linFun,x0,xdata,ydata1);

ydata2 = log(data(1:5,3));
y_hat2 = lsqcurvefit(@linFun,x0,xdata,ydata2);

ydata3 = log(data(1:5,4));
y_hat3 = lsqcurvefit(@linFun,x0,xdata,ydata3);

x = 1:1000;
m = y_hat1(1); b = y_hat1(2);
figure
hold on;
plot(x,m*x + b)
plot(xdata,ydata1,'o')

m = y_hat2(1); b = y_hat2(2);
plot(x,m*x + b)
plot(xdata,ydata2,'o')
legend('S Normal fit','S Normal data', 'S Tumor fit','S Tumor data')
title('Linear fitting on logarithm of Data')

% 3) Use Levenberg-Marquardt Method
% used as described in the MATLAB documentation of lsqcurvefit
% function. We use the monoexponential function stated in the question
% description (F = expFun(x,data)).
x0 = [0 0];
ub = [14000 0.1];
lb = [1000 0];
xdata  = data(1:5,1);
ydata1 = data(1:5,2);
y_hat1 = lsqcurvefit(@expFun,x0,xdata,ydata1,lb,ub, ...
    'Algorithm = levenberg-marquardt')

x = 1:1000;
y_hat1P = y_hat1(1)*exp(-x.*y_hat1(2));

ydata2 = data(1:5,3);
y_hat2 = lsqcurvefit(@expFun,x0,xdata,ydata2,lb,ub, ...
    'Algorithm = levenberg-marquardt')
y_hat2P = y_hat2(1)*exp(-x.*y_hat2(2));

figure
hold on;
plot(x,y_hat1P)
plot(xdata,ydata1,'o')
plot(x,y_hat2P)
plot(xdata,ydata2,'o')
legend('Estimate - S Normal','Data - S Normal','Estimate - S Tumour','Data - S Tumour')
% 4) Compare plots on logarithmic y scale
x = 1:1000;
y_hat1P = y_hat1(1)*exp(-x.*y_hat1(2));
figure
semilogy(x,y_hat1P)
hold on;
semilogy(xdata,ydata1,'o')
semilogy(x,y_hat2P)
semilogy(xdata,ydata2,'o')
title('Logarithmic plot of the Levenberg-Marquardt Method')
legend('Estimate - S Normal','Data - S Normal','Estimate - S Tumour','Data - S Tumour')

% 5) Using non-linear curve fitting over the entire data set
x0 = [0 0];
ub = [14000 0.1];
lb = [1000 0];
xdata  = data(:,1);
ydata1 = data(:,2);
ydata2 = data(:,3);

% Remove outliers [actually has little effect here]
xdata = [xdata(1:13); xdata(15:end)];
ydata1 = [ydata1(1:13); ydata1(15:end)];
ydata2 = [ydata2(1:13); ydata2(15:end)];

[y_hat1,resnorm1,residual1] = lsqcurvefit(@expFun,x0,xdata,ydata1,lb,ub, ...
    'Algorithm = levenberg-marquardt')
[y_hat2,resnorm2,residual2] = lsqcurvefit(@expFun,x0,xdata,ydata2,lb,ub, ...
    'Algorithm = levenberg-marquardt')

x = 1:3500;
y_hat1P = y_hat1(1)*exp(-x.*y_hat1(2));
y_hat2P = y_hat2(1)*exp(-x.*y_hat2(2));
figure
hold on;
plot(x,y_hat1P)
plot(xdata,ydata1,'o')
plot(x,y_hat2P)
plot(xdata,ydata2,'o')
title('Complete plot of the Levenberg-Marquardt Method')
legend('Estimate - S Normal','Data - S Normal','Estimate - S Tumour','Data - S Tumour')