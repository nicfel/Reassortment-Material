% builds a Coevo xml from the sequence data
clear; fclose('all');
% get all fasta files
segs_files = dir('data/*.fasta');

% define from when to when to take sequence
from = 1990;
to = 2020;

% set the random number generator
rng(48937345);

% define the number of sequences to use 
% nrsequences = 400;
nrsequences = 200;


seqs = cell(0,0);
seq_seqs = cell(0,0);
segments = cell(0,0);
for i = 1 : length(segs_files)
    fasta = fastaread(['data/' segs_files(i).name]);
    seq_seqs{i} = cell(0,0);
    tmp = strsplit(segs_files(i).name, '_');
    tmp = strsplit(tmp{2}, '.');
    segments{i} = tmp{1};
    for j = 1 : length(fasta)
        seqs{end+1} = fasta(j).Header;
        seq_seqs{i}{j} = fasta(j).Header;
    end
end

% get all unique sequences
unique_seqs = unique(seqs);


% delete any sequence for which the data isn't specified to the day
year = zeros(0,0);
for i = length(unique_seqs) : -1 : 1
    tmp = strsplit(unique_seqs{i}, '|');
    tmp = strsplit(tmp{2}, '-');
    year(i) = str2double(tmp{1});
    if length(tmp) < 3
        unique_seqs(i) = [];   
        year(i) = [];
    end        
end
% check if the dimension of the year vector is correct
if length(year)~=length(unique_seqs)
    error('error in the definition of the number of years');
end
% get the weights of each sequence as the inverse number of samples per
% year if the year is within from to
weights = zeros(length(year),1);
for i = 1 : length(weights)
    if year(i) >= from && year(i) <= to
        weights(i) = 1/(sum(year==year(i)));
    end
end

use_seqs = zeros(min(sum(weights>0), nrsequences),1);
for j = 1 : length(use_seqs)    
    add_seq = randsample(length(unique_seqs), 1, true, weights);
    year(add_seq) = -1;

    % recompute weights
    weights = zeros(length(year),1);
    for i = 1 : length(weights)
        if year(i) >= from && year(i) <= to
            weights(i) = 1/(sum(year==year(i)));
        end
    end
    use_seqs(j) = add_seq;
end


system('rm -r pairxmls');
system('mkdir pairxmls');

for a = 1 : length(segments)
    for b = a+1 : length(segments)
        use_segs{1} = segments{a};
        use_segs{2} = segments{b};
        seg_ind(1) = a;
        seg_ind(2) = b;

        for r = 0 : 2
            f = fopen('../template.xml');
            g = fopen(['pairxmls/h1n1pdm_' segments{a} '_' segments{b} '_rep' num2str(r) '.xml'], 'w');
            while ~feof(f)
                line = fgets(f);
                if ~isempty(strfind(line, 'insert_alignment'))
                    for i = 1 : length(segments)
                        fprintf(g, '\t<data id="%s">\n',segments{i});
                        fasta = fastaread(['data/' segs_files(i).name]);

                        seq_length(i) = length(fasta(1).Sequence);

                        for j = 1 : length(use_seqs)
                            ind = find(ismember(seq_seqs{i}, unique_seqs{use_seqs(j)}));
                            fprintf(g, '\t\t<sequence id="%s.%s" taxon="%s" totalcount="4" value="%s"/>\n',...
                                 segments{i}, unique_seqs{use_seqs(j)},...
                                 unique_seqs{use_seqs(j)}, strrep(fasta(ind).Sequence, 'U','T'));
                        end
                        fprintf(g, '\t</data>\n');

                    end
                elseif ~isempty(strfind(line, 'insert_run_header'))
                     fprintf(g, '\t\t<run id="mcmc" spec="beast.coupledMCMC.CoupledMCMC" logHeatedChains="true" chainLength="2500000" storeEvery="1000000" deltaTemperature="0.05" chains="4" resampleEvery="5000">\n');            
%                      fprintf(g, '\t\t<run id="mcmc" spec="beast.coupledMCMC.CoupledMCMC" logHeatedChains="true" chainLength="2500000" storeEvery="1000000" deltaTemperature="0.025" chains="4" resampleEvery="5000">\n');            
                elseif ~isempty(strfind(line, 'insert_taxa'))
                    for j = 1 : length(use_seqs)
                        fprintf(g, '\t\t\t<taxon spec="Taxon" id="%s"/>\n', unique_seqs{use_seqs(j)});
                    end
                elseif ~isempty(strfind(line, 'insert_sampling_times'))
                    for j = 1 : length(use_seqs)
                        time = strsplit(unique_seqs{use_seqs(j)}, '|');
                        tmp = strsplit(time{2}, '-');
                        % put the sampling time mid month 
                        if length(tmp) == 2
                            tmp{3} = '15';
                            time{2} = [time{2} '-15'];
                        end
                        deztime = (datenum(time{2},'yyyy-mm-dd')- datenum(tmp{1},'yyyy'))...
                            /(datenum(num2str(str2double(tmp{1})+1),'yyyy')-datenum(tmp{1},'yyyy'))...
                            +str2double(tmp{1});
                        dezstring = sprintf('%.8f', deztime);
                        fprintf(g, '\t\t\t\t%s=',unique_seqs{use_seqs(j)});
                        if j < length(use_seqs)
                            fprintf(g, '%s,\n', dezstring);
                        else
                            fprintf(g, '%s\n', dezstring);
                        end
                    end  
                elseif ~isempty(strfind(line, 'insert_nr_segments'))
                    fprintf(g, strrep(line, 'insert_nr_segments', num2str(length(use_segs))));

                elseif ~isempty(strfind(line, 'insert_parameters'))
                    disp(use_segs)
                    for s = 1 : length(use_segs)
                        fprintf(g, '\t\t\t\t\t\t<parameter id="kappa.s:%s_1" lower="0.0" name="stateNode">%f</parameter>\n',use_segs{s}, exprnd(1));
                        fprintf(g, '\t\t\t\t\t\t<parameter id="kappa.s:%s_3" lower="0.0" name="stateNode">%f</parameter>\n',use_segs{s}, exprnd(1));
                        fprintf(g, '\t\t\t\t\t\t<parameter id="mutationRate.s:%s_1" name="stateNode">1</parameter>\n',use_segs{s});
                        fprintf(g, '\t\t\t\t\t\t<parameter id="mutationRate.s:%s_3" name="stateNode">1</parameter>\n',use_segs{s});
                        fprintf(g, '\t\t\t\t\t\t<parameter id="gammaShape.s:%s_1" name="stateNode">%f</parameter>\n',use_segs{s}, exprnd(1));
                        fprintf(g, '\t\t\t\t\t\t<parameter id="gammaShape.s:%s_3" name="stateNode">%f</parameter>\n',use_segs{s}, exprnd(1));
                        fprintf(g, '\t\t\t\t\t\t<parameter id="freqParameter.s:%s_1" dimension="4" lower="0.0" name="stateNode" upper="1.0">0.25</parameter>\n',use_segs{s});
                        fprintf(g, '\t\t\t\t\t\t<parameter id="freqParameter.s:%s_3" dimension="4" lower="0.0" name="stateNode" upper="1.0">0.25</parameter>\n',use_segs{s});
                    end
                elseif ~isempty(strfind(line, 'insert_priors'))
                    for s = 1 : length(use_segs)

                        fprintf(g, '\t\t\t\t\t\t\t\t<prior id="KappaPrior.s:%s_1" name="distribution" x="@kappa.s:%s_1">\n', use_segs{s}, use_segs{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t<LogNormal id="LogNormalDistributionModel.%s.1_1" name="distr" M="1.0" S="1.25"/>\n', use_segs{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t</prior>\n');
                        fprintf(g, '\t\t\t\t\t\t\t\t<prior id="KappaPrior.s:%s_3" name="distribution" x="@kappa.s:%s_3">\n', use_segs{s}, use_segs{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t<LogNormal id="LogNormalDistributionModel.%s.1_3" name="distr" M="1.0" S="1.25"/>\n', use_segs{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t</prior>\n');
                        fprintf(g, '\t\t\t\t\t\t\t\t<prior id="GammaPrior.s:%s_1" name="distribution" x="@gammaShape.s:%s_1">\n', use_segs{s}, use_segs{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t<Exponential id="ExponentialDistribution.%s_3" name="distr"/>\n', use_segs{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t</prior>\n');
                        fprintf(g, '\t\t\t\t\t\t\t\t<prior id="GammaPrior.s:%s_3" name="distribution" x="@gammaShape.s:%s_3">\n', use_segs{s}, use_segs{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t<Exponential id="ExponentialDistribution.%s_30" name="distr"/>\n', use_segs{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t</prior>\n');
                    end
%                         fprintf(g, '\t\t\t\t\t\t\t\t<prior id="ReassortmentPrior" name="distribution" x="@reassortmentRate">\n');
%                         fprintf(g, '\t\t\t\t\t\t\t\t\t<LogNormal id="LogNormalDistributionModel.Reassortment" name="distr" M="0" S="4"/>\n');
%                         fprintf(g, '\t\t\t\t\t\t\t\t</prior>\n');
                elseif ~isempty(strfind(line, 'insert_operators'))
                    for s = 1 : length(use_segs)
                        fprintf(g, '\t\t\t\t<operator id="KappaScaler.s:%s_1" spec="ScaleOperator" parameter="@kappa.s:%s_1" scaleFactor="0.5" weight="0.1"/>\n', use_segs{s}, use_segs{s});
                        fprintf(g, '\t\t\t\t<operator id="KappaScaler.s:%s_3" spec="ScaleOperator" parameter="@kappa.s:%s_3" scaleFactor="0.5" weight="0.1"/>\n', use_segs{s}, use_segs{s});
                        fprintf(g, '\t\t\t\t<operator id="FrequenciesExchanger.s:%s_1" spec="DeltaExchangeOperator" delta="0.01" weight="0.1" parameter="@freqParameter.s:%s_1"/>\n', use_segs{s}, use_segs{s});
                        fprintf(g, '\t\t\t\t<operator id="FrequenciesExchanger.s:%s_3" spec="DeltaExchangeOperator" delta="0.01" weight="0.1" parameter="@freqParameter.s:%s_3"/>\n', use_segs{s}, use_segs{s});
                        fprintf(g, '\t\t\t\t<operator id="alpha_scaler_1.%s" spec="ScaleOperator" parameter="@gammaShape.s:%s_1" scaleFactor="0.75" weight="0.1"/>\n', use_segs{s}, use_segs{s});
                        fprintf(g, '\t\t\t\t<operator id="alpha_scaler_3.%s" spec="ScaleOperator" parameter="@gammaShape.s:%s_3" scaleFactor="0.75" weight="0.1"/>\n', use_segs{s}, use_segs{s});
                    end
                 elseif ~isempty(strfind(line, 'insert_mut_par'))
                    for s = 1 : length(use_segs)                   
                        fprintf(g, '\t\t\t\t\t<parameter idref="mutationRate.s:%s_1"/>\n', use_segs{s});
                        fprintf(g, '\t\t\t\t\t<parameter idref="mutationRate.s:%s_3"/>\n', use_segs{s});
                    end
                 elseif ~isempty(strfind(line, 'insert_param_log'))
                    for s = 1 : length(use_segs)                   
                        fprintf(g, '\t\t\t\t<log idref="kappa.s:%s_1"/>\n', use_segs{s});
                        fprintf(g, '\t\t\t\t<log idref="kappa.s:%s_3"/>\n', use_segs{s});
                        fprintf(g, '\t\t\t\t<log idref="mutationRate.s:%s_1"/>\n', use_segs{s});
                        fprintf(g, '\t\t\t\t<log idref="mutationRate.s:%s_3"/>\n', use_segs{s});
                        fprintf(g, '\t\t\t\t<log idref="gammaShape.s:%s_1"/>\n', use_segs{s});
                        fprintf(g, '\t\t\t\t<log idref="gammaShape.s:%s_3"/>\n', use_segs{s});
                        fprintf(g, '\t\t\t\t<log idref="freqParameter.s:%s_1"/>\n', use_segs{s});
                        fprintf(g, '\t\t\t\t<log idref="freqParameter.s:%s_3"/>\n', use_segs{s});
                    end
                elseif ~isempty(strfind(line, 'insert_tree_likelihood'))
                    for i = 1 : length(use_segs)
                        fprintf(g, '\t\t\t\t\t\t\t<distribution id="treeLikelihood.%s_1" spec="ThreadedTreeLikelihood" tree="@%s.tree" useAmbiguities="true">\n',use_segs{i},use_segs{i});
                        fprintf(g, '\t\t\t\t\t\t\t\t<data id="%s_1" spec="FilteredAlignment" filter="1::3,2::3">\n',use_segs{i});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t<data idref="%s"/>\n',use_segs{i});
                        fprintf(g, '\t\t\t\t\t\t\t\t</data>\n');
                        fprintf(g, '\t\t\t\t\t\t\t\t<siteModel id="SiteModel.s:%s_1" spec="SiteModel"  gammaCategoryCount="4" shape="@gammaShape.s:%s_1" mutationRate="@mutationRate.s:%s_1">\n',use_segs{i},use_segs{i},use_segs{i});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t<substModel id="hky.s:%s_1" spec="HKY" kappa="@kappa.s:%s_1">\n',use_segs{i},use_segs{i});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t\t<frequencies id="estimatedFreqs.s:%s_1" spec="Frequencies" frequencies="@freqParameter.s:%s_1"/>\n',use_segs{i},use_segs{i});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t</substModel>\n');
                        fprintf(g, '\t\t\t\t\t\t\t\t</siteModel>\n');
                        fprintf(g, '\t\t\t\t\t\t\t\t<branchRateModel id="StrictClock.%s_1" spec="beast.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c"/>\n',use_segs{i});
                        fprintf(g, '\t\t\t\t\t\t\t</distribution>\n');
                        fprintf(g, '\t\t\t\t\t\t\t<distribution id="treeLikelihood.%s_3" spec="ThreadedTreeLikelihood" tree="@%s.tree" useAmbiguities="true">\n',use_segs{i},use_segs{i});
                        fprintf(g, '\t\t\t\t\t\t\t\t<data id="%s_3" spec="FilteredAlignment" filter="3::3">\n',use_segs{i});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t<data idref="%s"/>\n',use_segs{i});
                        fprintf(g, '\t\t\t\t\t\t\t\t</data>\n');
                        fprintf(g, '\t\t\t\t\t\t\t\t<siteModel id="SiteModel.s:%s_3" spec="SiteModel"  gammaCategoryCount="4" shape="@gammaShape.s:%s_3" mutationRate="@mutationRate.s:%s_3">\n',use_segs{i},use_segs{i},use_segs{i});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t<substModel id="hky.s:%s_3" spec="HKY" kappa="@kappa.s:%s_3">\n',use_segs{i},use_segs{i});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t\t<frequencies id="estimatedFreqs.s:%s_3" spec="Frequencies" frequencies="@freqParameter.s:%s_3"/>\n',use_segs{i},use_segs{i});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t</substModel>\n');
                        fprintf(g, '\t\t\t\t\t\t\t\t</siteModel>\n');
                        fprintf(g, '\t\t\t\t\t\t\t\t<branchRateModel id="StrictClock.%s_3" spec="beast.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c"/>\n',use_segs{i});
                        fprintf(g, '\t\t\t\t\t\t\t</distribution>\n');
                    end
                elseif ~isempty(strfind(line, 'insert_weights'))
                    for i = 1 : length(use_segs)
                        fprintf(g, '%d %d ', round(seq_length(i)*2/3),round(seq_length(i)/3));
                    end
                elseif ~isempty(strfind(line, 'insert_seg_tree_state'))
                    for i = 1 : length(use_segs)
                         fprintf(g, '\t\t\t\t\t\t<stateNode id="%s.tree" spec="Tree" taxonset="@taxonSet" trait="@traitSet"/>\n', use_segs{i});                 
                    end
                elseif ~isempty(strfind(line, 'insert_seg_tree'))
                    for i = 1 : length(use_segs)
                         fprintf(g, '\t\t\t\t\t<segmentTree idref="%s.tree"/>\n', use_segs{i});                 
                    end
                elseif ~isempty(strfind(line, 'insert_seg_logger'))
                    for i = 1 : length(use_segs)
                        fprintf(g, '\t\t\t<logger spec="Logger" logEvery="50000" mode="tree" fileName="$(filebase).%s.trees">\n', use_segs{i});           
                        fprintf(g, '\t\t\t\t<log idref="%s.tree"/>\n', use_segs{i});           
                        fprintf(g, '\t\t\t</logger>\n');
                    end

                elseif ~isempty(strfind(line, 'insert_stats_log'))
                    for i = 1 : length(use_segs)
                        fprintf(g, '\t\t\t\t<log spec="TreeStatLogger" tree="@%s.tree"/>\n', use_segs{i});      
                    end
                elseif ~isempty(strfind(line, 'insert_seg_tree'))
                    for i = 1 : length(use_segs)
                         fprintf(g, '\t\t\t\t\t<segmentTree idref="%s.tree"/>\n', use_segs{i});                 
                    end

                else
                    fprintf(g, line);
                end
            end
            fclose(f);fclose(g);
        end   
    end
end
