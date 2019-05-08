% analyse reassortment distance distributions
clear
f = fopen('combined/infB.distance.txt');
reassortmentNode = cell(9,9);
reassortmentMRCA = cell(9,9);

count = 1;
while ~feof(f)
    line = strsplit(fgets(f));
    disp(count)
    
    if length(line)>0
        % get the mean distance of reassortment
        for i = 1 : 8
            for j = i+1 :9
                reassortmentNode{i,j}{count} = zeros(0,0);
                reassortmentMRCA{i,j}{count} = zeros(0,0);
            end
        end
        for i = 2 : length(line)-1
            tmp = strsplit(line{i}, ',');
            for j = 1 : length(tmp)
                tmp2 = strsplit(tmp{j}, ':');
                tmp3 = strsplit(tmp2{1},'-');
                segs = sort(str2double(tmp3))+1;
                reassortmentNode{segs(1),segs(2)}{count}(end+1,1) = str2double(tmp2{2});
                reassortmentMRCA{segs(1),segs(2)}{count}(end+1,1) = str2double(tmp2{3});
            end
        end
    end    
    count = count+1;
end



%% 
x = [5:0.01:20];

segments = {'HA', 'MP', 'NA', 'NP', 'NS', 'PA', 'PB1', 'PB2', 'prior'};

threshold_mrca = 30;
threshold_node = 20;

c = 1;
for i = 1 %: 8
    for j = 1:8
        if j>i
            y = zeros(0,1);
            y_prior = zeros(0,0);
            for k = 1 : length(reassortmentMRCA{i,j})
                indices = reassortmentMRCA{i,j}{k}>threshold_mrca;
                y = [y ; min(reassortmentNode{i,j}{k}(indices))];
                indices_prior = reassortmentMRCA{i,9}{k}>threshold_mrca;
                y_prior = [y_prior ;min(reassortmentNode{i,9}{k}(indices_prior))];
            end
            
            
                       

            subplot(1,8,c)
            [yval,xval] = ksdensity(y,x); 
            [yval_prior,xval_prior] = ksdensity(y_prior,x); 
            plot(xval_prior,yval_prior, 'red');
            hold on
            plot(xval,yval, 'blue'); hold on
%             alpha(.5)
            hold off
            title(segments{j})
            ylabel(segments{i})
%             xlabel('log(sum RD pair/sum RD empty)')
        end
        c = c+1;
    end
end
