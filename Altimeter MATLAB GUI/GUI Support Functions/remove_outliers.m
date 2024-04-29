function [T, P] = remove_outliers(t,p,sig)
P = [];
T = [];
ii = 1;
mean_p = mean(p);
standard_dev = std(p);
for i = 1:length(p)
    if (abs(p(i)-mean_p) < sig*standard_dev) && (p(i) > 0)
        P(ii) = p(i);
        T(ii) = t(i);
        ii = ii +1;
    end
end
end