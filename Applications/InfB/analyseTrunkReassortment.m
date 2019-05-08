% analyse if reassortment events are more or less likely to occur on the
% trunk
clear
f = fopen(['combined/infB.trunk.txt']);

count = 1;
probs = zeros(0,2);
while ~feof(f)
    line = str2double(strsplit(strtrim(fgets(f))));

    if length(line)>0
        probs(count,1) = log((line(1)/line(3))/(line(2)/line(4)));
%             probs(count,2) = line(2)/line(4);
    end

    count = count+1;
end
ksdensity(probs(:,1));hold on
%     ksdensity(probs(:,2));hold on
%     xlim([0, ])
xlabel('prob of reassortment on trunk over prob of reassortment off trunk');

