function [prob] = getBinomPdf(data, p)
if p<=0 || p>=1
    prob = -9999999999;
    return;
end
prob = 0;
for i = 1 : length(data.before)
    add_prob = log(binopdf(data.after(i), data.before(i), p) + ...
    	binopdf(data.before(i)-data.after(i), data.before(i), p));
    prob = prob + add_prob;
end
end
