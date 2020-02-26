function [ f, sigma, cb, diff, rsq, sse ] = errorfunctionFATNAV( data )
% Fits errorfunction to input data. The fit, value of sigma, confidence
% bounds of sigma (cb), and 2 measures for the goodness of fit (rsquare,
% sse) can be returned as output.

x = double(data(:,1));
y = double(data(:,2));
sp = [min(y), max(y)-min(y), 1, 0];
% fo = fitoptions('Method', 'NonlinearLeastSquares','StartPoint', sp,...
%     'Algorithm', 'Levenberg-Marquardt');

fo = fitoptions('Method', 'NonlinearLeastSquares','StartPoint', sp,...
    'Algorithm', 'Levenberg-Marquardt', 'Robust','Bisquare');
ft = fittype(@(l,h,sig,mu,x) l+h*0.5*erf((x-mu)./(sqrt(2)*sig)), 'options',fo);
[f, gof]  = fit(x,y,ft);

param = coeffvalues(f);
sigma = param(3);
conf = confint(f);
cb = conf(:,3);
rsq = gof.rsquare;
sse = gof.sse;
diff = param(2);

end

