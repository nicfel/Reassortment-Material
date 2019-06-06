clear

currdir = cd;

virus = {'h1n1pdm',...
        'h1n1sea',...
        'h3n2',...
        'infB',...
        'h3n2new'};
workingdir = {'H1N1pandemic', ...
        'H1N1seasonal', ...
        'H3N2', ...
        'InfB',...
        'H3N2'};

from = [2009, ...
        1990, ...
        1980,...
        1920,...
        2015];
to =   [2019, ...
        2010, ...
        2020, ...
        2019,...
        2020];


temperature = [0.005, 0.005, 0.005, 0.005, 0.005];
temperature = temperature*2;

% system('rm -r xmls');
% system('mkdir xmls');

% values for the random number generator that were at some point sampled at
% random from randi(1000000,100,1). This is done to be able to reproduce the
% exact subsampled datasets
rng_values = [738601,620765,812970,106938,585114];

max_year_offset = [1,1,1,1,1];

for i = 5 : length(virus)
    disp(virus{i})
    getXMLallsegmentsTipEstimation(virus{i}, workingdir{i}, 20, from(i), to(i), temperature(i), rng_values(i), 100, max_year_offset(i));
    cd(currdir)
    getJointCoalallsegmentsTipEstimate(virus{i}, workingdir{i}, 20, from(i), to(i), temperature(i), rng_values(i), 100, max_year_offset(i));
    cd(currdir)
end