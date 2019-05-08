% analyse if reassortment events are more or less likely to occur on the
% trunk
clear
name = dir('combined/*trunk.txt');
mapping = [1,2,5,3,4];
for i = 1 : length(name)
    f = fopen(['combined/' name(i).name]);

    count = 1;
    probs = zeros(0,1);
    while ~feof(f)
        line = str2double(strsplit(strtrim(fgets(f))));

        if length(line)>0
            probs(count,1) = (line(1)/line(3))-(line(2)/line(4));
%             probs(count,2) = line(2)/line(4);
        end

        count = count+1;
    end
    disp(sum(probs>0.0)/length(probs))
    subplot(length(name),1,mapping(i))
    ksdensity(probs(:,1));hold on
    xlim([-0.5, 0.5])
    title(name(i).name)
end

count = 1;
values = ce
probs = zeros(0,1);
for i = 2 : length(name)
    f = fopen(['combined/' name(i).name]);

    while ~feof(f)
        line = str2double(strsplit(strtrim(fgets(f))));
        if length(line)>0
            probs(count,1) = (line(1)/line(3))-(line(2)/line(4));
        end
        count = count+1;

    end
end


disp(sum(probs>0.0)/length(probs))
subplot(length(name),1,mapping(i))
ksdensity(probs(:,1));hold on
xlim([-0.5, 0.5])
title(name(i).name)




% name = {'','_old','_recent'};
% titlename = {'1980-2010','1980-1995','1995-2010'};
% for i = 1 : 3
%     f = fopen(['simulation/trunk' name{i} '_reassortment.txt']);
% 
%     count = 1;
%     probs = zeros(0,1);
%     while ~feof(f)
%         line = str2double(strsplit(strtrim(fgets(f))));
% 
%         if length(line)>0
%             probs(count,1) = log((line(1)/line(3))/(line(2)/line(4)));
% %             probs(count,2) = line(2)/line(4);
%         end
% 
%         count = count+1;
%     end
%     
%     probs(isinf(probs)) = 0;
%     
%     disp(sum(probs>0.0)/length(probs))
%     subplot(3,1,i)
%     ksdensity(probs(:,1));hold on
% %     ksdensity(probs(:,2));hold on
%     xlim([-1, 1])
%     title(titlename{i})
%     legend('prob rea on trunk over prob rea off trunk');
% end
% 
% 
