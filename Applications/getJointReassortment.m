% compute probabilities of splitting/joinig segments
clear
dirs = {'H5N1', 'H3N2', 'InfB'};
time = zeros(0,0);
run_nr=0;
skip_nr = 200;

joint = cell(0,0);
split = cell(0,0);

nr_before = cell(0,0);
nr_after = cell(0,0);
for r = 1:length(dirs)
    c = 1;
    skip_count = 1;


    logfiles = dir([dirs{r} '/out/*.reassortment.log']);
    f = fopen([dirs{r} '/out/' logfiles(1).name]);fgets(f);
    while skip_count>skip_nr
        fgets(f);
        skip_nr = skip_nr+1;
    end
    while ~feof(f)
%         disp(run_nr)
        run_nr = run_nr+1;
        line = fgets(f);
        splitline = strsplit(line, '\t');
        joint{c} = zeros(8,8);
        split{c} = zeros(8,8);

        for rea = 2 : length(splitline)
            tmp = strsplit(splitline{rea},':');
            vals = strrep(tmp{1}, '{','');
            vals = strrep(vals, '}','');
            
            % get the segments before a reassortment event
            before = str2double(strsplit(vals, ','));
            
            vals = strrep(tmp{2}, '{','');
            vals = strrep(vals, '}','');

            % get the segments after a reassortment event
            after1 = str2double(strsplit(vals, ','));
            
            vals = strrep(tmp{3}, '{','');
            vals = strrep(vals, '}','');
            
            % get the segments after a reassortment event
            after2 = str2double(strsplit(vals, ','));

                        
            for a = 1 : 7
                for b = a+1 : 8
                    %% check if both are on the same lineage before
                    if sum(ismember(before,a-1))==1 && sum(ismember(before,b-1))==1
                        if sum(ismember(after1,a-1))==1 && sum(ismember(after1,b-1))==1 || ...
                            sum(ismember(after2,a-1))==1 && sum(ismember(after2,b-1))==1
                            joint{c}(a,b) = joint{c}(a,b) + 1;
                        else
                            split{c}(a,b) = split{c}(a,b) + 1;
                        end
                    end
                end
            end
        end
        c=c+1;

    end
    fclose(f);
    
    figure()
    joint_vals = zeros(size(split,2),8,8);
    split_vals = zeros(size(split,2),8,8);
    for i = 1 : size(split,2)
        if ~isempty(joint{i})
            joint_vals(i,:,:) = joint{i};
            split_vals(i,:,:) = split{i};
        else
            break;
        end
    end

    segments = {'HA', 'M1', 'NA', 'NP', 'NS1', 'PA', 'PB1', 'PB2'};
    c = 1;
    for a = 1:8
        for b = 1:8        
            if b>a
                subplot(8,8,c)
                [y,x]=ksdensity(joint_vals(:,a,b));
                [y2,x2]=ksdensity(split_vals(:,a,b));
                plot(x,y, 'Color', 'red');hold on
                plot(x2,y2, 'Color', 'blue');
                title(sprintf('%s %s', segments{a}, segments{b}))
% %                 xlim([-2 2])
            end
            c = c+1;
        end
    end
end
