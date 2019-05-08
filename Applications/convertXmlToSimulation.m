function [] = convertXmlToSimulation(folder)
% converts an xml to run a reassorty analyses to perform simulations using
% the same tip dates and the mean estimated reassortment rate and effective
% population size
    cd(folder)
    log_files = dir('combined/*.log');
    for i = 1 : length(log_files)
        t = importdata(['combined/' log_files(i).name]);
        % look for the reassormtment and pop size headers
        rea_ind = -1; pop_ind = -1;
        for j = 1 : length(t.textdata)
            if strcmp(t.textdata{j}, 'reassortmentRate')
                rea_ind = j;
            elseif strcmp(t.textdata{j}, 'popSize.t')
                pop_ind = j;
            end
        end

        mean_rea = mean(t.data(:,rea_ind));
        mean_pop = mean(t.data(:,pop_ind));

        % get the taxa and dates from the xmls corresponding to the log files
        f = fopen(['xmls/' strrep(log_files(i).name, '.log', '.xml')]);
        taxa = cell(0,0);
        dates = cell(0,0);
        while ~feof(f)
            line = fgets(f);
            if contains(line, '<taxonset id="taxonSet" spec="TaxonSet">')
                line = fgets(f);
                while ~contains(line, '</taxonset>')
                    taxa{end+1} = line;
                    line = fgets(f);
                end           
            elseif contains(line, '<trait spec="TraitSet" traitname="date-forward" id="traitSet">')
                line = fgets(f);
                while ~contains(line, '<taxa idref="taxonSet"/>')
                    dates{end+1,1} = line;
                    line = fgets(f);
                end
            end
        end
        fclose(f);
        % build a simulation file with that info
        f = fopen('../simulation_template.xml');
        g = fopen(['simulation/' strrep(log_files(i).name, 'rep0.log', 'sim.xml')], 'w');
        while ~feof(f)
            line = fgets(f);
            if contains(line, 'insert_Ne')
                fprintf(g, strrep(line, 'insert_Ne', num2str(mean_pop)));
            elseif contains(line, 'insert_reassortment_rate')
                fprintf(g, strrep(line, 'insert_reassortment_rate', num2str(mean_rea)));
            elseif contains(line, 'insert_taxa')
                for j = 1 : length(taxa)
                    fprintf(g, '%s', taxa{j});
                end
            elseif contains(line, 'insert_dates')
                for j = 1 : length(taxa)
                    fprintf(g, '%s', dates{j});
                end
            else
                fprintf(g, line);
            end        
        end
        fclose('all');
    end
end