function [e,v] = myltnormal(ltnmean,ltnvar,ltnpoint)
if ltnvar <= 10^-30
    ltnvar = 10^-30;
end
ltnstd = sqrt(ltnvar);
stdpoint = (ltnpoint - ltnmean)/ltnstd;
% cdf_stdpoint = 1-cdf('norm',stdpoint,0,1);
% pdf_stdpoint = pdf('norm',stdpoint,0,1);
% 
% 
% if (cdf_stdpoint==0)
%     e1 = 0;
%     v1 = 0;
% else
%     p = pdf_stdpoint/cdf_stdpoint;
%     e1 = ltnmean+p*ltnstd;
%     v1 = (1+stdpoint*p-p^2)*ltnvar;
% end

t = erfc(stdpoint/sqrt(2));
if (t == 0)
    e = 0;
    v = 0;
else
    p = exp(0.5*log(2/pi)-0.5*stdpoint^2-log(t));
    e = ltnmean+p*ltnstd;
    v = (1+stdpoint*p-p^2)*ltnvar;
end