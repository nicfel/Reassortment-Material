% compare the pairwise phylogenetic distances of simulated vs. inferred
clear;fclose('all');
% get all the simulated files
simtrees = dir('simout/*Seg1.tree');
% define the number of segements
nrsegments = 4;
% file to print the results to
f_false = fopen('noevent_detection.csv','w');
fprintf(f_false,'max,difference,detection_probability,run\n');

true_points = zeros(0,3);
for i = 1:length(simtrees)
    disp(simtrees(i).name)
    % read in the simulated segments
    for j = 1 : nrsegments
        [sim.pdist{j}, sim.names{j}] = ...
            getPairwiseDistances(['simout/' strrep(simtrees(i).name,...
            'Seg1',sprintf('Seg%d', j))],1);
        namebase = strrep(simtrees(i).name,'Seg1',sprintf('Seg%d', j));
        namebase = strrep(namebase, 'trueS','s');
        [sim.nucdiff{j}] = ...
            getNucleotideDifferences(['simout/' strrep(namebase, '.tree','.alignment.nexus')]);

    end
    % read in the inferred segment trees
    inf_name = strrep(simtrees(i).name, 'trueS','s');
    inf_name = strrep(inf_name, 'sim_','inf_');
    inf_name = strrep(inf_name, '.tree','.combined.trees');
    for j = 1 : nrsegments        
        [inf.pdist{j}, inf.names{j}] = ...
            getPairwiseDistances(['combined/' strrep(inf_name,...
            'seg1',sprintf('seg%d', j))],1);        
    end
    
    %%
    % find leafs between which there was a reassortment event
    cutoff = 0.000001;
    has_reassortment = false(size(sim.pdist{1}{1}));
    max_dist = zeros(size(sim.pdist{1}{1}));
    max_nucdiff = zeros(size(sim.pdist{1}{1}));
    min_nucdiff = zeros(size(sim.pdist{1}{1}));
    min_nucdiff = zeros(size(sim.pdist{1}{1}));
    max_pairwise_dist = zeros(size(sim.pdist{1}{1}));

    for a = 1 : size(sim.pdist{1}{1},1)
        for b = a+1 : size(sim.pdist{1}{1},2)
            vals = zeros(nrsegments,1);
            nucdiff = zeros(nrsegments,1);
            for j = 1:nrsegments
                vals(j) = sim.pdist{j}{1}(a,b);
                nucdiff(j) = sim.nucdiff{j}(a,b);
            end
            if max(abs(diff(vals)))>cutoff
                has_reassortment(a,b) = 1;
                max_dist(a,b) = floor(max(abs(diff(vals)))*10^6)/10^6;
           end
            max_nucdiff(a,b) = max(nucdiff);
            min_nucdiff(a,b) = min(nucdiff);
       end
    end
    
    
    
    
    % see if the inference detects a reassortment event seperating the two
    reassortment_prob = zeros(size(sim.pdist{1}{1}));
    seperation = zeros(size(sim.pdist{1}{1}));
    for a = 1 : size(sim.pdist{1}{1},1)
        for b = a+1 : size(sim.pdist{1}{1},2)
            is_reassorted_count = 0;
            sep_val = zeros(0,0);
            for j = 1 : length(inf.pdist{1})
                for k = 1:nrsegments
                    vals(k) = inf.pdist{k}{j}(a,b);
                end
                if max(abs(diff(vals)))>cutoff
                    is_reassorted_count = is_reassorted_count+1;
                    sep_val(end+1) = max(abs(diff(vals)));
                end
            end
            reassortment_prob(a,b) = is_reassorted_count/length(inf.pdist{1});  
            if length(sep_val)>0
                mean_seperation(a,b) = mean(sep_val);
            else
                mean_seperation(a,b) = 0;
            end
        end
    end
    % get the run number
    runval = strsplit(simtrees(i).name,'_');
    runval = strsplit(runval{2},'.');
    runval = runval{1};
    
    
        % try to get the unique events
    uni_max_dist = unique(max_dist);
    
    for a = 2 : length(uni_max_dist)
        indices = find(max_dist==uni_max_dist(a));
        % find the indice with the lowest pairwise distance
        [~,ind_min] = min(max_nucdiff(indices));
        true_points(end+1,:) = [max_nucdiff(indices(ind_min)),min_nucdiff(indices(ind_min)), reassortment_prob(indices(ind_min))];
    end

    
    % print to file
    for a = 1 : size(sim.pdist{1}{1},1)
        for b = a+1 : size(sim.pdist{1}{1},2)
            if has_reassortment(a,b)
%                 fprintf(f_true, '%f,%f,%f,%s\n',max_nucdiff(a,b),min_nucdiff(a,b), reassortment_prob(a,b), runval);
            else
                fprintf(f_false, '%f,%f,%f,%s\n',max_nucdiff(a,b),min_nucdiff(a,b), reassortment_prob(a,b),runval);
            end
       end
    end
end
%% if there is more than one point with the same coordinates, print the mean
f_true = fopen('event_detection.csv','w');
fprintf(f_true,'max,min,detection_probability,run\n');

uni_x = unique(true_points(:,1));
uni_y = unique(true_points(:,2));
for a = 1 : length(uni_x)
    for b = 1 : length(uni_y)
        indices = find(true_points(:,1)==uni_x(a) & true_points(:,2)==uni_y(b));
        if length(indices)>0
            fprintf(f_true, '%f,%f,%f,%s\n',uni_x(a),uni_y(b), mean(true_points(indices,3)), runval);        
           
        end
    end
end
fclose('all');