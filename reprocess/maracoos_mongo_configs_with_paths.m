% Rutgers HFRadar Processing Toolbox
%
% reprocess_with_mongo_configs.m
%
% This script connects to the MongoDB database and loads the config file
% described by the 'config_name' variable. 
%
% Created by Mike Smith (michaesm@marine.rutgers.edu) on 4/13/2019
addpath(genpath('.'))
init_time = now;

% Start and end times of range you want to reprocess
start_time = datetime(2017,12,1,0,0,0);
end_time = datetime(2017,12,1,2,0,0);

config_name = 'maracoos_long_range_reprocess'; % Name of configs document
collection = 'totals_configs'; % MongoDB collection
conn = mongo_database; % Connect to MongoDB database


t0_str = string(start_time, 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z''');
t1_str = string(end_time, 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z''');

% query document. document contains configs in the format required by hfrprogs
document = find(conn, collection, 'Query', ['{"name": "'  config_name '"}']);
conf = document.conf;
conf.Radials.RangeBearSlop = repmat(1e-10, [numel(conf.Radials.Sites), 2]);

site_strings = join(compose('{"Site": "%s"}', cell2mat(conf.Radials.Sites)), ', ');
query_site_str = sprintf('{$and: [{$or: [%s]}, ', site_strings{1});
query_time_str = sprintf('{"TimeStampStr": {$gt: "%s"}}, {"TimeStampStr": {$lte: "%s"}}]}', t0_str, t1_str);
query_str = [query_site_str, query_time_str];

radial_table = find(conn, 'radials','Query', query_str, 'Projection', '{"_id":0, "Site": 1, "Path": 1}');
radial_table = struct2table(radial_table);

% Build hourly timestamps
time_steps = end_time:-1/24:start_time;  

for x = 1:1:length(time_steps)
    t = time_steps(x);
    t_str = datestr(t, 'yyyy_mm_dd_HH00.ruv');

    % create subset of table that only includes time t
%     radial_table_sub = radial_table(strcmp(t, radial_table.TimeStampStr),:);
    radial_table_sub = radial_table(contains(radial_table.Path, t_str),:);
    
    % Process current timestamp
    fprintf(1, '****************************************\n');
    fprintf(1, '  Current time: %s\n',datestr(now));
    fprintf(1, '  Processing data time: %s\n', t_str);

    % Hourly Total Creation
    try
        fprintf(1, 'Starting Totals_driver\n');
        codar_driver_totals(datenum(t), conf, radial_table_sub.Path);
    catch
        fprintf(1, 'Driver_totals failed because: \n');
        res = lasterror;                
        fprintf(1, '%s\n',res.message);
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