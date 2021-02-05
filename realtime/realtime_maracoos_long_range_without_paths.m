% Rutgers HFRadar Processing Toolbox
%
% realtime_maracoos_long_range.m
%
% This is the main input script for total generation. This script works
% for all monostatic SeaSonde system types. Bistatic will be implemented in
% the future.
%
% This script connects to a MongoDB database to gather a json document that
% contains configuration data in the format required for HFRProgs. From the
% configuration document, the script grabs a table of all radial files 
% during the past n hours (defined in the config document). The script then
% loops through each timestamp during the last n number of hours checking 
% if the config file wants to process Ideal, Masured, or Best Choice
% (operators choice of best pattern), or all four. A list of radial strings
% is passed to the create_hourly_totals function where the script checks
% if the data needs to be reprocessed (due to new radials)
% 
%
% Created by Mike Smith (michaesm@marine.rutgers.edu) on 5/21/2019
% See also radials2totals, CODAR_configuration, CODAR_driver_totals
clear, close all
init_time = now;

fprintf('CheckTotals.m Start Time: %s\n', datestr(init_time)); 
addpath(genpath('/home/codaradm/operational_scripts_test/totals_toolbox/totals'))
addpath(genpath('/home/codaradm/HFR_Progs-2_1_3beta'))

which maketotalsOI
which tuvOI

% Name of the MongoDB document containing configs you want to use
config_name = 'realtime_long_range';
infer = false; % infer radials to use based on available files in directory. 
              % assumes user has appropriate radial types for each site 
              
% |Name of the collection in the MongoDB that contains the processing
collection = 'totals_configs';

% Connect to MongoDB databas
conn = mongo_database;

% Find the configs document
c = find(conn, 'totals_configs', 'Query', ['{"name": "'  config_name '"}']);
c.conf.Radials.RangeBearSlop = repmat(1e-10, [numel(c.conf.Radials.Sites), 2]);

t1 = dateshift(datetime('now', 'TimeZone', 'UTC'), 'end', 'hour');
t0 = t1 - c.conf.Radials.HoursAgoToCheck/24;

%t1 = datetime(2020, 3, 18, 11, 0, 0)
%t0 = datetime(2020, 3, 11, 18, 0, 0)


% Build hourly timestamps
time_steps = datenum(t1):-1/24:datenum(t0);  
sites = c.conf.Radials.Sites;

for x = 1:1:length(time_steps)
      
    t = time_steps(x);

    % Process current timestamp
    fprintf(1, '****************************************\n');
    fprintf(1, '  Current time: %s\n', datestr(now));
    fprintf(1, '  Processing data time: %s\n', t);
    fprintf(1,  '****************************************\n');
    
    % Build cell arrays of potential ideal and measured radial files
    ideal_files = filenames_standard_filesystem(c.conf.Radials.BaseDir, c.conf.Radials.Sites(:), 'RDLi', t, c.conf.Radials.MonthFlag, c.conf.Radials.TypeFlag);
    measured_files = filenames_standard_filesystem(c.conf.Radials.BaseDir, c.conf.Radials.Sites(:), 'RDLm', t, c.conf.Radials.MonthFlag, c.conf.Radials.TypeFlag);
    
    if infer == true
        infer_files = [ideal_files; measured_files];
        existing = {};
        for f = 1:length(infer_files)
            if exist(infer_files{f}, 'file')
                existing = [existing; infer_files{f}];
            end
        end
        create_hourly_totals(t, c, existing, 'inferred');
    else
        if c.conf.Radials.ProcessByPattern.Ideal
            create_hourly_totals(t, c, ideal_files, 'Ideal', conn);
        end

        if c.conf.Radials.ProcessByPattern.Measured
            create_hourly_totals(t, c, measured_files, 'Measured', conn);
        end

        if c.conf.Radials.ProcessByPattern.BestChoice
            % Build cell array of potential best pattern selected radial files
            best_files = filenames_standard_filesystem(c.conf.Radials.BaseDir, c.conf.Radials.Sites(:), c.conf.Radials.Types(:), t, c.conf.Radials.MonthFlag, c.conf.Radials.TypeFlag);

            create_hourly_totals(t, c, best_files, 'BestChoice', conn);
        end
    end
    
end

fprintf(1, 'Total Creation Finished.\n');
fprintf(1, '---------------------------------------------------------------------------\n');

% Display script ending time and elapsed time to file.
fprintf(1, 'reprocess.m End Time: %s \n', datestr(now));
end_time = abs(etime(floor(datevec(init_time)), floor(datevec(now))));
elapsed_minutes = floor(end_time/60);
elapsed_seconds = mod(end_time, 60);
fprintf(1, 'Elapsed time is %s minutes and %s seconds.\n', num2str(elapsed_minutes), num2str(elapsed_seconds));
