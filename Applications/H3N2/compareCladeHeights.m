% compares heights of clades between the coalescent with and without
% reassortment
clear


logs = dir('trees/clades/*new*.trees');

g = fopen('trees/heights.csv', 'w');        
fprintf(g, 'post1,post2,lower1,lower2,upper1,upper2,median1,median2,segment,time\n');

for i = 1 : 1
    for j = 1 : length(logs)
        segment = strrep(logs(j).name, 'h3n2new.','');
        segment = strrep(segment, '.trees','');
        
        if i==1
            filename = strrep(logs(j).name, 'new','');
            time='all';
        else
            filename = logs(j).name;
            time='new';
        end
        
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
                % print the posterior values
                fprintf(g, '%s,%s,', strrep(post1, ']', ''), strrep(post2, ']', ''));
                
                % get the height hpd's
                fprintf(g, '%s,', strrep(tmp1{8}, 'height_95%_HPD={', ''));
                fprintf(g, '%s,', strrep(tmp2{8}, 'height_95%_HPD={', ''));
                
                fprintf(g, '%s,', strrep(tmp1{9}, '}', ''));
                fprintf(g, '%s,', strrep(tmp2{9}, '}', ''));   
                
                fprintf(g, '%s,', strrep(tmp1{10}, 'height_median=', ''));
                fprintf(g, '%s,', strrep(tmp2{10}, 'height_median=', ''));

                
                fprintf(g, '%s,%s\n', segment, time);
                
            else
                if contains(meta2{i}, 'posterior=')
                    error('nodes don''t match');
                end
            end
        end
    end
end
fclose('all')