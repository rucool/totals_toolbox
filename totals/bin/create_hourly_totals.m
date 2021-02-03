function create_hourly_totals(t, configs, file_list, pattern_field, conn)
%     tnum = datenum(t, 'yyyy-mm-ddTHH:MM:SSZ');
    tstr = datestr(t, 'xyyyymmddTHHMMSS');
    
    available_files_list = cell(0);
    for y = 1:size(file_list,1)
       if exist(char(file_list(y)), 'file')
           available_files_list  = [available_files_list; file_list(y) ];
       end
    end
    
    if ~exist('conn', 'var')
        % conn parameter does not exist... default to false
        codar_driver_totals(datenum(t), configs.conf, file_list, pattern_field);
    else
        if ~isfield(configs.totals_generated.(pattern_field), tstr)
            % Record does not exist. Create Totals and insert record into db.
            fprintf(1, 'A record for %s does not exist. Total Creation Starting. \n', t);
            old_radials = 0;
        else % Record exists, but check if the new number of radials is greater than the old number of radials.
            fprintf(1, 'A record for %s exists. Determining if re-processing is needed. \n', t);
            old_radials = configs.totals_generated.(pattern_field).(tstr).num_radials;
        end

        new_radials = length(available_files_list); % Number of Radials that are available for a specific time period

        if new_radials > old_radials
            % Hourly Total Creation
            try
                fprintf(1,  'Making Hourly Totals. \n');
                codar_driver_totals(datenum(t), configs.conf, file_list, pattern_field);

                if ~isempty(conn)
                    field = sprintf('totals_generated.%s.%s', pattern_field, datestr(t, 30));
                    
                    join_files = join(compose('"%s"', string(file_list)), ', ');
                    update_str = sprintf('[%s]', join_files{1});

                    % Add mongo insert total stuff below this line
                    find_query = sprintf('{"name":"%s"}', configs.name);
                    update_query = sprintf('{$set: {"%s.num_radials": %d, "%s.files": %s}}', field, new_radials, field, update_str);
                    update(conn, 'totals_configs', find_query, update_query);
                end
            catch
                fprintf(1, 'Total Creation failed because: \n');
                res = lasterror;
                fprintf(1, '%s\n',res.message);
            end    
        elseif new_radials < 2
            fprintf(1, 'Not enough radials to generate totals. Continuing to next timestep. \n');
        else
            fprintf(1, 'Record exists and timestamp does NOT need to be re-processed. \n');
        end
    end
end