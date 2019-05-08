% analyse reassortment distance distributions
clear
f = fopen('combined/h3n2.distance.txt');
reassortmentDist = cell(9,9);
count = 1;
while ~feof(f)
    line = strsplit(fgets(f));
    
    if length(line)>0
        % get the mean distance of reassortment
        for i = 1 : 8
            for j = i+1 :9
                reassortmentDist{i,j}{count} = zeros(0,0);
            end
        end
        for i = 2 : length(line)-1
            tmp = strsplit(line{i}, ',');
%             val = 0;
            for j = 1 : length(tmp)
                tmp2 = strsplit(tmp{j}, ':');
%                 val = val + str2double(tmp2{2});
                tmp3 = strsplit(tmp2{1},'-');
                segs = sort(str2double(tmp3))+1;
                reassortmentDist{segs(1),segs(2)}{count}(end+1,1) = str2double(tmp2{3})-str2double(tmp2{2});
            end
        end
    end    
    count = count+1;
end



%% 
x = [-20:1:20];

segments = {'HA', 'MP', 'NA', 'NP', 'NS', 'PA', 'PB1', 'PB2', 'prior'};

threshold=2;

c = 1;
for i = 1 : 8
    for j = 1:8
        if j>i
            y = zeros(0,1);
            y_prior = zeros(0,0);
            for k = 1 : length(reassortmentDist{i,j})
                y = [y ; sum(reassortmentDist{i,j}{k}>threshold)-sum(reassortmentDist{i,9}{k}>threshold)];
                y_prior = [y_prior ;sum(reassortmentDist{i,9}{k}>threshold)];
            end
            
            
                       

            subplot(7,8,c)
            [yval,xval] = hist(y,x); 
            [yval_prior,xval_prior] = hist(y_prior,x); 
            bar(xval,yval, 'blue'); hold on
            plot([0,0],[0 50], 'red', 'LineWidth',3)
%             bar(xval_prior,yval_prior, 'red')
%             alpha(0.5)
%             title(segments{j})
%             ylabel(segments{i})
%             xlabel('log(sum RD pair/sum RD empty)')
        end
        c = c+1;
    end
end
