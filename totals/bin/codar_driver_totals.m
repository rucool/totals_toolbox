function codar_driver_totals_qc(dtime, conf, F, pattern_field, varargin)
%
% Rutgers HFRadar Processing Toolbox
%
% Usage: CODAR_driver_totals(dtime, conf, varargin)
%
% This script processes radial files for input into total generation. 
%
% First, it converts the conf variable called from CODAR_configuration.m
% and converts it into a format that is useable by the HFRProgs toolbox.
% Then, it loads the radials from all available sites. After loading, it
% cleans the radials of bad data and masks the radials that are lying over
% the land. After this process, it calls makeTotalsOI.m to create the total
% vectors. Once the totals are generated, this script also cleans and masks
% out any bad data that may have been generated during total creation.
% 
% Results are written to the data directory in matlab formats
%
% See also CheckTotals, CODAR_configuration, CODAR_driver_totals,
% makeTotalsOI

% Convert conf variable into array that HFRProfs can handle
conf = HFRPdriver_default_conf( conf );

% Mandatory parameters required by user
mand_params = {'Radials.Sites', 'Radials.Types', 'Radials.RangeLims',...
               'Radials.BearLims', 'Totals.DomainName', 'Totals.GridFile'};
            
% Check if mandatory input arguments are satisfied
conf = checkParamValInputArgs( conf, {}, mand_params, varargin{:} );

% Load Radial Data
fprintf(1, 'Loading Radial Data. \n');
fprintf(1, '--------------------------------------------------------- \n');

Rorig = loadRDLFile_all(F(:));
fprintf(1, '--------------------------------------------------------- \n');

% If files contain no data, filter them out. 
ii = false(size(Rorig));

for j = 1:numel(Rorig)
    ii(j) = numel(Rorig(j).U) == 0;
end

missingRadials.FileNames = [ Rorig(ii).FileName ];
[missingRadials.TimeStamps,missingRadials.Sites,missingRadials.Types]...
      = parseRDLFileName( missingRadials.FileNames );

Rorig(ii) = [];
if isempty(Rorig)
    error( 'No data at this timestep.' );
end

% Get rid of stuff for missing files - if maskfile, rangelims or bearlims
% are missing, just don't do anything
try conf.Radials.MaskFiles(ii) = []; end
try conf.Radials.RangeLims(ii,:) = []; end
try conf.Radials.BearLims(ii,:) = []; end

if conf.QC.Process  % If conf.QC.Process is set to true
    % Remove radials that were flagged with a QC tests
    Rqc = cleanQCedRadials(Rorig, 'all: QC*', 4);
    
    % Remove radials that were flagged with a QC tests (header metadata)
    % Rqch = cleanQCedRadialFiles(Rqc, 'all: qc_qartod_*', 7);
    Rqch = cleanQCedRadialFiles(Rqc, 'QC09', 4);
    
    % Remove radials above MaxRadSpeed threshold
    fprintf(1, 'Cleaning Radials. \n');
    Rclean = cleanRadials( Rqch, conf.Radials.MaxRadSpeed );
else % Else conf.QC.Process is false
    fprintf(1, 'Cleaning Radials. \n');
    Rclean = cleanRadials(Rorig, conf.Radials.MaxRadSpeed );
end
    
% Remove radials over land
fprintf('Masking Radials. \n');
Rmask = maskRadials( Rclean, conf.Radials.MaskDir, 0);

%-------------------------------------------------
% Interpolate missing radials
fprintf(1, 'Interpolating Radials. \n');
fprintf(1, '--------------------------------------------------------- \n');
for n = 1:numel(Rmask)
    % Allow for possibilty of RangeLims and/or BearLims to be not defined.
    try
        RL = conf.Radials.RangeLims(n,:);
    catch
        RL = [];
    end
    try
        BL = conf.Radials.BearLims(n,:);
    catch
        BL = [];
    end
        
    % If there is only one radial, or maybe some other conditons that I'm 
    % not thinking of, then the interpolation will fail.  Set up a
    % try/catch block to keep things from failing.
    try
        Rinterp(n) = interpRadials(Rmask(n), ...
                            'RangeLims', RL, ...
                            'BearLims', BL, ...
                            'RangeDelta', conf.Radials.RangeBearSlop(n,1),...
                            'BearDelta', conf.Radials.RangeBearSlop(n,2),...
                            'MaxRangeGap', conf.Radials.RangeGap,...
                            'MaxBearGap', conf.Radials.BearGap,...
                            'CombineMethod', 'average');
    catch
        fprintf(1, 'Warning: ## interpRadials failed for Site: %s\n', Rmask(n).SiteName);
        res=lasterror;
        fprintf(1, '%s\n',res.message)
        Rinterp(n) = Rmask(n);
        Rinterp(n).ProcessingSteps{end+1} = 'Interpolation failed, revert to uninterpolated';
    end
    
    % Check for case of interpolation creating all NaN's.  Use 90 % as the 
    % threshold.  Replace with uninterpolated radials and warn the user.
    if (sum(~isnan(Rinterp(n).U)) < sum(~isnan(Rmask(n).U)) * 0.9)
        fprintf(1, '%s:\n',char(Rinterp(n).FileName));
        fprintf(1, 'probably not interpolated properly ... using uninterpolated data instead\n')
        tmp = Rinterp(n).ProcessingSteps;
        Rinterp(n) = Rmask(n);
        Rinterp(n).ProcessingSteps = tmp;
        Rinterp(n).ProcessingSteps{end+1} = 'Revert to uninterpolated';
    end
end
fprintf(1, '--------------------------------------------------------- \n');

% Load the total grid
[grid, fn, c] = loadDataFileWithChecks(conf.Totals.GridFile);
if c >= 100
  error( 'Could not find totals grid.' );
end

%% Make Unweighted Least Squares Totals
if conf.Totals.UWLS.Process
    fprintf(1, 'Generating UWLS totals. \n');
    [TUVuwls, RTUV] = makeTotals(Rinterp,...
        'Grid', grid,...
        'TimeStamp', dtime,...
        'spatthresh', conf.Totals.UWLS.spatthresh,...
        'tempthresh', 1/24/2-eps,...
        'DomainName', conf.Totals.DomainName,...
        'CreationInfo', conf.Totals.CreationInfo);
    
    TUV = TUVuwls;
    TUVqc = TUVuwls;
    
    
    [TUVclean, I] = cleanTotals(TUV,...
        conf.Totals.MaxTotSpeed);
    fprintf('%d totals removed by cleanTotals\n', sum(I(:)>0))

    % cleanTotals sets I=1 where speed is > maxspeed and I=0 for good values
    % so have to adjust to match QARTOD convention
    I(I == 1) = 4;
    I(I == 0) = 1;
    TUVqc.QC16 = I;
    clear I;
    
    % Mask totals    
    [~, IL]=maskTotals(TUV,conf.Totals.MaskFile,false);
    fprintf(1, '%d totals masked out. \n',sum(~IL(:)));
    
    % IL = logical index the same size as the original number of grid points
    % with true for every grid point that was kept (i.e. not trimmed)
    % adjust to match QARTOD convention
    %IL(IL == 0) = 4; % this statement does NOT work because IL is a logical variable
    
    I = ones(size(IL,1),1);
    I(IL == 0) = 4;
    TUVqc.QC18 = I;  % there is no QARTOD location test but will call this QC18 for now
    clear I, IL;
    
    [~,I] = cleanTotals(TUVclean,conf.Totals.MaxTotSpeed,conf.Totals.cleanTotalsVarargin{:}); 
    % I = Index of data points that were removed.  This will have a zero for
    %     good data, a 1 for data whose speed is above maxspd, a 2 for data
    %     violated first error condition, a 4 for data violating second error
    %     condition, etc.  For data violating multiple conditions, the result
    %     will be a sum of all the violated conditions. 
    % b/c max speed already called index variable will not be a sum since both tests cannot fail
    I(I == 2) = 4;
    I(I == 0) = 1;
    TUVqc.QC15 = I;
    clear I;
    
    
    TUVqc.PRIM = max([TUVqc.QC15, TUVqc.QC16,TUVqc.QC18 ],[],2);
    minqc = min([TUVqc.QC15, TUVqc.QC16,TUVqc.QC18 ],[],2);
    minqc1 = find(minqc == 1);
    prim2 = find(TUVqc.PRIM == 2);
    setto1 = intersect(minqc1, prim2); 
    if ~isempty(setto1)
      TUVqc.PRIM(setto1) = 1;  % if at least one other test has passed, set the primary flag to pass instead of not evaluated
    end

    TUVqc.qc_operator_mask = ones(length(TUV.U),1).*2;

    % Save results
    [tdn, tfn] = datenum_to_directory_filename([conf.Totals.UWLS.BaseDir lower(pattern_field)],...
        dtime,...
        conf.Totals.UWLS.FilePrefix,...
        conf.Totals.UWLS.FileSuffix,...
        conf.MonthFlag);

    tdn = tdn{1};

    if ~exist(tdn, 'dir')
      mkdir(tdn);
    end

    save(fullfile(tdn, tfn{1}),'conf','missingRadials','RTUV','TUVuwls','TUV','TUVqc');
end
fprintf(1, '--------------------------------------------------------- \n');

if conf.Totals.OI.Process
    % Create OI Totals
    % Start with Rmask (pre interpolation)
    RTUV = Rmask;

    % Call makeTotalsOI to generate the totals.
    fprintf(1, 'Generating OI totals. \n');
    [TUVorig, RTUV]=makeTotalsOI(RTUV,...
        'Grid', grid,...
        'TimeStamp', dtime,...
        'mdlvar', conf.Totals.OI.mdlvar,...
        'errvar', conf.Totals.OI.errvar,...
        'sx', conf.Totals.OI.sx,...
        'sy', conf.Totals.OI.sy,...
        'tempthresh',conf.Totals.OI.tempthresh, ...
        'DomainName',conf.Totals.DomainName, ...
        'CreationInfo',conf.Totals.CreationInfo);
    
    TUV = TUVorig;
    TUVqc = TUVorig;
    
    % Clean totals
    [TUVclean, I] = cleanTotals(TUV, conf.Totals.MaxTotSpeed); %, ...

    fprintf(1, '%d totals removed by cleanTotals. \n',sum(I(:)>0));

    % cleanTotals sets I=1 where speed is > maxspeed and I=0 for good values
    % so have to adjust to match QARTOD convention
    I(I == 1) = 4;
    I(I == 0) = 1;
    TUVqc.QC16 = I;
    clear I;

    % Mask totals
    [~, IL]=maskTotals(TUV,conf.Totals.MaskFile,false);
    fprintf(1, '%d totals masked out. \n',sum(~IL(:)));
    
    % IL = logical index the same size as the original number of grid points
    % with true for every grid point that was kept (i.e. not trimmed)
    % adjust to match QARTOD convention
    %IL(IL == 0) = 4; % this statement does NOT work because IL is a logical variable
    
    I = ones(size(IL,1),1);
    I(IL == 0) = 4;
    TUVqc.QC18 = I;  % there is no QARTOD location test but will call this QC18 for now
    clear I, IL;
    
    [~,I] = cleanTotals(TUVclean,conf.Totals.MaxTotSpeed,conf.Totals.OI.cleanTotalsVarargin{1});  
    % I = Index of data points that were removed.  This will have a zero for
    %     good data, a 1 for data whose speed is above maxspd, a 2 for data
    %     violated first error condition, a 4 for data violating second error
    %     condition, etc.  For data violating multiple conditions, the result
    %     will be a sum of all the violated conditions.
    % b/c max speed already called index variable will not be a sum since both tests cannot fail
    I(I == 2) = 4;
    I(I == 0) = 1;
    TUVqc.QC19 = I;  % there is no QARTOD OI Uerr uncertainty test but will call this QC19 for now
    clear I;
    
    [~,I] = cleanTotals(TUVclean,conf.Totals.MaxTotSpeed,conf.Totals.OI.cleanTotalsVarargin{2});  
    % b/c max speed already called index variable will not be a sum since both tests cannot fail
    I(I == 2) = 4;
    I(I == 0) = 1;
    TUVqc.QC20 = I; % there is no QARTOD OI Uerr uncertainty test but will call this QC20 for now
    clear I;
    
    TUVqc.PRIM = max([TUVqc.QC16,TUVqc.QC18,TUVqc.QC19,TUVqc.QC20 ],[],2);
    minqc = min([TUVqc.QC16,TUVqc.QC18,TUVqc.QC19,TUVqc.QC20 ],[],2);
    minqc1 = find(minqc == 1);
    prim2 = find(TUVqc.PRIM == 2);
    setto1 = intersect(minqc1, prim2); 
    if ~isempty(setto1)
      TUVqc.PRIM(setto1) = 1;  % if at least one other test has passed, set the primary flag to pass instead of not evaluated
    end
  
    TUVqc.qc_operator_mask = ones(length(TUV.U),1).*2;
    
    % %% ----------------------------------------------------------------------
    % %% This section of code was a contribution from Erick Fredj to gap fill the data
    % % Masking is a destructive process, any totals current inside the mask will be
    % % set to zero.
    % 
    % fprintf(1, 'Starting Smooth Total Field. \n');
    % fprintf(1, '--------------------------------------------------------- \n');
    % 
    % mask=load(conf.OSN.BestCoverageFile);
    % hfrcvrg = inpolygon(TUV.LonLat(:,1),TUV.LonLat(:,2),mask(:,1),mask(:,2));
    % 
    % % Robust Smooth
    % TUVosn = TUV;
    % TUVosn.CreationInfo= 'Erick Fredj';
    % 
    % U=TUVosn.U(hfrcvrg);
    % V=TUVosn.V(hfrcvrg);
    % 
    % % set to reset TUVs.U to NaN
    % TUVosn.U = NaN(size(TUVosn.U));
    % % set to reset TUVs.V to NaN
    % TUVosn.V = NaN(size(TUVosn.V));
    % 
    % %% this function smoothn is located in toolbox_eric_fredj
    % Vs = smoothn({U,V},'robust');
    % 
    % TUVosn.U(hfrcvrg)=Vs{1};
    % TUVosn.V(hfrcvrg)=Vs{2};

    %%-------------------------------------------------------------------------

    % Save results
    [tdn, tfn] = datenum_to_directory_filename([conf.Totals.OI.BaseDir lower(pattern_field)],...
        dtime,...
        conf.Totals.OI.FilePrefix,...
        conf.Totals.FileSuffix,...
        conf.Totals.MonthFlag);
    
    tdn = tdn{1};

    if ~exist(tdn, 'dir')
      mkdir(tdn);
    end

    save(fullfile(tdn, tfn{1}),'conf','missingRadials','RTUV','TUVorig','TUV','TUVqc');
end
fprintf(1, '--------------------------------------------------------- \n');