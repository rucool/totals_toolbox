% Rutgers HFRadar Processing Toolbox
%
% CheckTotals.m
%
% This script connects to the MongoDB database to see if the system type
% (long, mid, standard) exists. If the system type exists, this script
% calls radials2totals.m, which by default, determines the times that 
% totals need to be generated. 
%
% Created by Mike Smith (michaesm@marine.rutgers.edu) on 4/13/2019
% See also codar_driver_totals
addpath(genpath('.'))
init_time = now;

% Start and end times of range you want to reprocess
start_time = datetime(2017,1,1,10,0,0);
end_time = datetime(2017,1,1,12,0,0);

% Name of the MongoDB document containing configs you want to use
config_name = 'maracoos_long_range_reprocess';
infer = true; % infer radials to use based on available files in directory. 
              % assumes user has appropriate radial types for each site 

% Connect to MongoDB database
conn = mongo_database;

% Find the configs document
c = find(conn, 'totals_configs', 'Query', ['{"name": "'  config_name '"}']);
c.conf.Radials.RangeBearSlop = repmat(1e-10, [numel(c.conf.Radials.Sites), 2]);

% Build hourly timestamps
time_steps = end_time:-1/24:start_time;  
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