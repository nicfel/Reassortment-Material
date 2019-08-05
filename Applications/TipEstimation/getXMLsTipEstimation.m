clear

currdir = cd;

virus = {'h3n2',...
        'h3n2new'};
workingdir = {'H3N2', ...
        'H3N2'};

from = [1995,...
        1995];
to =   [2019,...
        2019];
interval_length = [9, 1];
    


temperature = [0.005, 0.005, 0.005, 0.005, 0.005];
temperature = temperature*2;

% system('rm -r xmls');
% system('mkdir xmls');

% values for the random number generator that were at some point sampled at
% random from randi(1000000,100,1). This is done to be able to reproduce the
% exact subsampled datasets
rng_values = [738601,620765,812970,106938,585114];

% max_year_offset = [1,1,1,1,0];

for i = 1 : length(virus)
    disp(virus{i})
    getXMLallsegmentsTipEstimation(virus{i}, workingdir{i}, 20, from(i), to(i), temperature(i), rng_values(i), 500, 0, interval_length(i));
    cd(currdir)
    getJointCoalallsegmentsTipEstimate(virus{i}, workingdir{i}, 20, from(i), to(i), temperature(i), rng_values(i), 500, 0, interval_length(i));
    cd(currdir)
end
