% compute probabilities of splitting/joinig segments
clear
dirs = {'H5N1', 'H3N2', 'InfB'};
time = zeros(0,0);
run_nr=0;
skip_nr = 20;
skip_count = 1;
c = 1;

joint = zeros(8,8);
split = zeros(8,8);

for r = 1:length(dirs)
    nr_before = cell(0,0);
    nr_after = cell(0,0);
    c = 1;
    skip_count = 1;


    logfiles = dir([dirs{r} '/out/*.reassortment.log']);
    f = fopen([dirs{r} '/out/' logfiles(1).name]);fgets(f);
    while skip_count>skip_nr
        fgets(f);
        skip_nr = skip_nr+1;
    end
    while ~feof(f)
        disp(run_nr)
        run_nr = run_nr+1;
        line = fgets(f);
        splitline = strsplit(line, '\t');
        for rea = 2 : length(splitline)
            tmp = strsplit(splitline{rea},':');
            vals = strrep(tmp{1}, '{','');
            vals = strrep(vals, '}','');
            
            % get the segments before a reassortment event
            before = str2double(strsplit(vals, ','));
            
            vals = strrep(tmp{2}, '{','');
            vals = strrep(vals, '}','');

            % get the segments after a reassortment event
            after = str2double(strsplit(vals, ','));
            
            data(run_nr).before(rea-1) = length(before);
            data(run_nr).after(rea-1) = length(after);
            
%             for a = 1 : 7
%                 for b = a+1 : 8
%                     if sum(ismember(before,a-1))==1 && sum(ismember(before,b-1))==1
%                         if sum(ismember(after,a-1))==1 && sum(ismember(after,b-1))==1
%                             joint(a,b) = joint(a,b) + 1;
%                         else
%                             split(a,b) = split(a,b) + 1;
%                         end
%                     end
%                 end
%             end

            time(c) = str2double(tmp{end});
            c=c+1;
        end
%         break;
    end
    fclose(f);
end

