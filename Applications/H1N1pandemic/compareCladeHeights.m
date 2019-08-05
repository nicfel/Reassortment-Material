% compares heights of clades between the coalescent with and without
% reassortment
clear


logs = dir('trees/clades/*.trees');

g = fopen('trees/heights.csv', 'w');        
fprintf(g, 'post1,post2,lower1,lower2,upper1,upper2,median1,median2,segment,virus\n');

for i = 1 : 1
    for j = 1 : length(logs)
        tmp = strsplit(logs(j).name, '.');
        virus = tmp{1};
        segment = tmp{2};
        
        filename = logs(j).name;
        
        % get the first tree
        f = fopen(['trees/clades/' filename]);
        while ~feof(f)
            line = fgets(f);
            if length(line)>1000
                tree1 = line;
            end
        end
        
        % get the second tree
        fclose(f);f = fopen(['trees/mcc/' filename]);
        while ~feof(f)
            line = fgets(f);
            if length(line)>1000
                tree2 = line;
            end
        end
        fclose(f);
        % get all meta data values
        meta1 = regexp(tree1, '\[(.*?)\]', 'match');
        meta2 = regexp(tree2, '\[(.*?)\]', 'match');
        
        
        for k = 1 : length(meta1)
            if contains(meta1{k}, 'posterior=')
                % split the strings
                tmp1 = strsplit(meta1{k}, ',');
                tmp2 = strsplit(meta2{k}, ',');
                
                % sanity check
                if length(tmp1) ~= length(tmp2)
                    error('node meta data don''t match');
                end
                % get the posterior values
                post1 = strrep(tmp1{end}, 'posterior=', '');
                post2 = strrep(tmp2{end}, 'posterior=', '');
                
                % check if a posterior is 0
                if ~strcmp(post1,'0.0]') && ~strcmp(post2,'0.0]')

                    % print the posterior values
                    fprintf(g, '%s,%s,', strrep(post1, ']', ''), strrep(post2, ']', ''));

                    % get the height hpd's
                    fprintf(g, '%s,', strrep(tmp1{8}, 'height_95%_HPD={', ''));
                    fprintf(g, '%s,', strrep(tmp2{8}, 'height_95%_HPD={', ''));

                    fprintf(g, '%s,', strrep(tmp1{9}, '}', ''));
                    fprintf(g, '%s,', strrep(tmp2{9}, '}', ''));   

                    fprintf(g, '%s,', strrep(tmp1{10}, 'height_median=', ''));
                    fprintf(g, '%s,', strrep(tmp2{10}, 'height_median=', ''));


                    fprintf(g, '%s,%s\n', segment, virus);
                end
                
            else
                if contains(meta2{i}, 'posterior=')
                    error('nodes don''t match');
                end
            end
        end
    end
end
fclose('all');


%% compute the clade support for the mcc trees using both methods
g = fopen('trees/clades.csv', 'w');        
fprintf(g, 'post,segment,method,virus\n');

for j = 1 : length(logs)
    tmp = strsplit(logs(j).name, '.');
    virus = tmp{1};
    segment = tmp{2};
    
    filename = logs(j).name;
        
    % get the first tree
    f = fopen(['trees/mcc/' filename]);
    while ~feof(f)
        line = fgets(f);
        if length(line)>1000
            tree1 = line;
        end
    end

    % get the second tree
    fclose(f);f = fopen(['trees/mcc/' strrep(filename, virus,[virus 'norea'])]);
    while ~feof(f)
        line = fgets(f);
        if length(line)>1000
            tree2 = line;
        end
    end
    fclose(f);
    
    % get all meta data values
    meta = regexp(tree1, '\[(.*?)\]', 'match');

    for k = 1 : length(meta)
        if contains(meta{k}, 'posterior=')
            % split the strings
            tmp = strsplit(meta{k}, ',');

            % get the posterior values
            post1 = strrep(tmp{end}, 'posterior=', '');

            % print the posterior values
            fprintf(g, '%s,%s,rea,%s\n', strrep(post1, ']', ''), segment, virus);
        end
    end    
    % get all meta data values
    meta = regexp(tree2, '\[(.*?)\]', 'match');

    for k = 1 : length(meta)
        if contains(meta{k}, 'posterior=')
            % split the strings
            tmp = strsplit(meta{k}, ',');

            % get the posterior values
            post1 = strrep(tmp{end}, 'posterior=', '');

            % print the posterior values
            fprintf(g, '%s,%s,norea,%s\n', strrep(post1, ']', ''), segment, virus);
        end
    end    

end
fclose('all');



