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
start_time = now;

fprintf('CheckTotals.m Start Time: %s\n', datestr(start_time)); 
addpath(genpath('.'))

% Name of the MongoDB document containing configs you want to use
config_name = 'swarm_mid_range_realtime';

% |Name of the collection in the MongoDB that contains the processing
collection = 'totals_configs';

% Connect to MongoDB databas
conn = mongo_database;

% Get the document containing processing configs via name of the document
configs = find(conn, collection, 'Query', ['{"name": "'  config_name '"}']);
configs.conf.Radials.RangeBearSlop = repmat(1e-10, [numel(configs.conf.Radials.Sites), 2]);

% Get latest radial time and build past n number of hours
latest = find(conn, 'radials', 'Sort', '{"TimeStamp":-1}', 'Limit', 1, 'Projection','{"TimeStamp":1.0}');
t1 = datetime(latest.TimeStamp.x_date, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z''');
t0 = t1 - configs.conf.Radials.HoursAgoToCheck/24;
t0_str = string(t0, 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z''');
t1_str = string(t1, 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z''');

site_strings = join(compose('{"Site": "%s"}', cell2mat(configs.conf.Radials.Sites)), ', ');
query_site_str = sprintf('{$and: [{$or: [%s]}, ', site_strings{1});
query_time_str = sprintf('{"TimeStampStr": {$gt: "%s"}}, {"TimeStampStr": {$lte: "%s"}}]}', t0_str, t1_str);
query_str = [query_site_str, query_time_str];

radial_table = find(conn, 'radials','Query', query_str, 'Projection', '{ "PatternType": 1, "Site": 1, "Path": 1, "_id": 0, "TimeStampStr": 1}');
radial_table = struct2table(radial_table);
times = table2array(sortrows(unique(radial_table(:,4)),1, 'descend'));

% Close MongoDB Connection
close(conn)

for x = 1:1:length(times)
    % Connect to MongoDB database
    conn = mongo_database;
    t = times{x};
    % create subset of table that only includes time t
    radial_table_sub = radial_table(strcmp(t, radial_table.TimeStampStr),:);

    % Process current timestamp
    fprintf(1, '****************************************\n');
    fprintf(1, '  Current time: %s\n', datestr(now));
    fprintf(1, '  Processing data time: %s\n', t);
    fprintf(1,  '************e****************************\n');

    if configs.conf.Radials.ProcessByPattern.Ideal
        % list all ideal files from radial_table_sub
        ideal_files=radial_table_sub.Path(strcmp('ideal', radial_table_sub.PatternType));
        create_hourly_totals(datenum(t, 'yyyy-mm-ddTHH:MM:SSZ'), configs, ideal_files, 'Ideal', conn);
    end

    if configs.conf.Radials.ProcessByPattern.Measured
        % list all measured files from radial_table_sub
        measured_files=radial_table_sub.Path(strcmp('measured', radial_table_sub.PatternType));
        create_hourly_totals(datenum(t, 'yyyy-mm-ddTHH:MM:SSZ'), configs, measured_files, 'Measured', conn);
    end

    if configs.conf.Radials.ProcessByPattern.BestChoice
        % determine which sites are measured and which are ideal for "best" option
        measured_sites=configs.conf.Radials.Sites(strcmp('RDLm', configs.conf.Radials.Types));
        ideal_sites=configs.conf.Radials.Sites(strcmp('RDLi', configs.conf.Radials.Types));

        % list best files (ideal)
        best_files=radial_table_sub.Path(strcmp('ideal',radial_table_sub.PatternType)&ismember(radial_table_sub.Site,ideal_sites));

        % add best files (measured) to best files with ideal
        best_files=[best_files;radial_table_sub.Path(strcmp('measured',radial_table_sub.PatternType)&ismember(radial_table_sub.Site,measured_sites))];
        create_hourly_totals(datenum(t, 'yyyy-mm-ddTHH:MM:SSZ'), configs, best_files, 'BestChoice', conn);
    end
    
    % Close MongoDB Connection
    close(conn);
end


fprintf(1, 'Total Creation Finished.\n');
fprintf(1, '---------------------------------------------------------------------------\n');

% Display script ending time and elapsed time to file.
end_time = now;
fprintf(1, 'check_totals_long_range.m End Time: %s \n', datestr(end_time));
elapsed = abs(etime(datevec(end_time), datevec(start_time)));
elapsed_minutes = floor(elapsed/60);
elapsed_seconds = mod(elapsed, 60);
fprintf(1, 'Elapsed time is %s minutes and %s seconds.\n', num2str(elapsed_minutes), num2str(elapsed_seconds));
