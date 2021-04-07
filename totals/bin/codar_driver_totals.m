function codar_driver_totals(dtime, conf, F, pattern_field, varargin)
%
% Rutgers HFRadar Processing Toolbox
%
% Usage: CODAR_driver_totals(dtime, conf, varargin)
%
% This script processes radial files for input into total generation. 
%
% First, it converts the conf variable sourced from total_configs collection in the Mongo database
% and converts it into a format that is useable by the HFRProgs toolbox.
% Then, it loads the radials from all available sites. After loading, it
% cleans the radials of flagged data and masks the radials that are lying over
% the land. After this process, it calls makeTotalsOI.m to create the total
% vectors. Once the totals are generated, this script also flags any bad data
% that may have been generated during total creation.
% 
% Results are written to the data directory in matlab formats
%
% See also CheckTotals, CODAR_driver_totals, makeTotalsOI
%  and total_configs collection in the Mongo database

% Convert conf variable into array that HFRProfs can handle
conf = HFRPdriver_default_conf( conf );

% Mandatory parameters required by user
mand_params = {'Radials.Sites', 'Radials.Types', 'Radials.RangeLims',...
               'Radials.BearLims', 'Totals.DomainName', 'Totals.GridFile'};
            
% Check if mandatory input arguments are satisfied
conf = checkParamValInputArgs( conf, {}, mand_params, varargin{:} );

% do not make a folder called BestChoice, pattern_field will only be used to create a folder 
% if the choice is to compute totals with all ideal or all measured
if strcmp(pattern_field,'BestChoice')
   pattern_field = '';  
end

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
end
    
% Remove radials over land
% In case the QARTOD valid location test, didn't flag all of the radials over land, this
% mask will ensure those vectors are not used in the total vector calculation
fprintf('Masking Radials. \n');
Rmask = maskRadials( Rqc, conf.Radials.MaskDir, 0);

%-------------------------------------------------

% Load the total grid
[grid, fn, c] = loadDataFileWithChecks(conf.Totals.GridFile);
if c >= 100
  error( 'Could not find totals grid.' );
end

%% Make Unweighted Least Squares Totals
if conf.Totals.UWLS.Process
    fprintf(1, 'Generating UWLS totals. \n');
    [TUVuwls, RTUV] = makeTotals(Rmask,...
        'Grid', grid,...
        'TimeStamp', dtime,...
        'spatthresh', conf.Totals.UWLS.spatthresh,...
        'tempthresh', 1/24/2-eps,...
        'DomainName', conf.Totals.DomainName,...
        'CreationInfo', conf.Totals.CreationInfo);
    
    TUV = TUVuwls;
    TUVqc = TUVuwls;
    
    
    [TUVclean, I] = cleanTotals(TUV,...
        conf.QC.QC303_MaxTotSpeed);
    fprintf('%d totals removed by cleanTotals\n', sum(I(:)>0))

    % cleanTotals sets I=1 where speed is > maxspeed and I=0 for good values
    % so have to adjust to match QARTOD convention
    I(I == 1) = 4;
    I(I == 0) = 1;
    TUVqc.QC303 = I;
    clear I;
    
    % Mask totals    
    [~, IL]=maskTotals(TUV,conf.QC.QC305_MaskFile,false,0,0);
    fprintf(1, '%d totals masked out. \n',sum(~IL(:)));
    
    % IL = logical index the same size as the original number of grid points
    % with true for every grid point that was kept (i.e. not trimmed)
    % adjust to match QARTOD convention
    %IL(IL == 0) = 4; % this statement does NOT work because IL is a logical variable
    
    I = ones(size(IL,1),1);
    I(IL == 0) = 4;
    TUVqc.QC305 = I;  % there is no QARTOD location test but will call this QC305 for now
    clear I, IL;
    
    GDOP_test_parameters = {conf.QC.GDOP_UWLS.QC302_GDOP.ErrorTypeStr conf.QC.GDOP_UWLS.QC302_GDOP.ErrorVariableStr conf.QC.GDOP_UWLS.QC302_GDOP.MaxError}';

    [~,I] = cleanTotals(TUVclean,conf.QC.QC303_MaxTotSpeed,GDOP_test_parameters); 
    % I = Index of data points that were removed.  This will have a zero for
    %     good data, a 1 for data whose speed is above maxspd, a 2 for data
    %     violated first error condition, a 4 for data violating second error
    %     condition, etc.  For data violating multiple conditions, the result
    %     will be a sum of all the violated conditions. 
    % b/c max speed already called index variable will not be a sum since both tests cannot fail
    I(I == 2) = 4;
    I(I == 0) = 1;
    TUVqc.QC302 = I;
    clear I;
    
    TUVqc.qc_operator_flag = ones(length(TUV.U),1).*2;
    
    TUVqc.PRIM = max([TUVqc.QC302, TUVqc.QC303,TUVqc.QC305,TUVqc.qc_operator_flag],[],2);
    minqc = min([TUVqc.QC302, TUVqc.QC303,TUVqc.QC305, TUVqc.qc_operator_flag ],[],2);
    minqc1 = find(minqc == 1);
    prim2 = find(TUVqc.PRIM == 2);
    setto1 = intersect(minqc1, prim2); 
    if ~isempty(setto1)
      TUVqc.PRIM(setto1) = 1;  % if at least one other test has passed, set the primary flag to pass instead of not evaluated
    end

   
    TUV = TUVqc;
    TUV.MARACOOS_TUV_struct_version = 'undefined';

    TUVmetadata.conf = conf;
    TUVmetadata.radial_metadata = [];
    for rr = 1:size(RTUV,1)
      clear tmp;
      tmp = char(RTUV(rr).FileName);
      TUVmetadata.radial_metadata{rr} = tmp(max(strfind(tmp,'/'))+1:end);
    end
    TUVmetadata.radial_num_sites = size(RTUV,1);
    TUVmetadata.missingRadials = missingRadials;
    TUVmetadata.attributes = [];
    TUVmetadata.header = ['QCPrimaryFlagDefinition: Highest flag value of these test flags: QC303, QC305, QC302, operator (Note: Primary flag will be set to 2 only if ALL tests were not evaluated)'];


    % Save results
    [tdn, tfn] = datenum_to_directory_filename([conf.Totals.UWLS.BaseDir,...
        'mat/', conf.SystemType,'/', lower(pattern_field)],...
        dtime,...
        conf.Totals.UWLS.FilePrefix,...
        conf.Totals.UWLS.FileSuffix,...
        conf.MonthFlag);

    tdn = tdn{1};
 
    if ~exist(tdn, 'dir')
      mkdir(tdn);
    end

    save(fullfile(tdn, tfn{1}),'RTUV','TUV','TUVmetadata');
    clear tdn tfn
    
end
fprintf(1, '--------------------------------------------------------- \n');

if conf.Totals.OI.Process
    % Create OI Totals
    RTUV = Rmask;

    if strcmp(conf.Totals.OI.weighting_function, 'exponential')
      wf = 2;
    elseif strcmp(conf.Totals.OI.weighting_function, 'gaussian')
      wf = 1;
    else
      warning('Incorrect setting for choice of OI weighting function.  It should be set to exponential or gaussian.') 
    end
      
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
        'weighting', wf, ...
        'DomainName',conf.Totals.DomainName, ...
        'CreationInfo',conf.Totals.CreationInfo);
    
    TUV = TUVorig;
    TUVqc = TUVorig;
    
    % Clean totals
    [TUVclean, I] = cleanTotals(TUV, conf.QC.QC303_MaxTotSpeed); %, ...

    fprintf(1, '%d totals removed by cleanTotals. \n',sum(I(:)>0));

    % cleanTotals sets I=1 where speed is > maxspeed and I=0 for good values
    % so have to adjust to match QARTOD convention
    I(I == 1) = 4;
    I(I == 0) = 1;
    TUVqc.QC303 = I;
    clear I;

    % Mask totals
    [~, IL]=maskTotals(TUV,conf.QC.QC305_MaskFile,false,0,0);
    fprintf(1, '%d totals masked out. \n',sum(~IL(:)));
    
    % IL = logical index the same size as the original number of grid points
    % with true for every grid point that was kept (i.e. not trimmed)
    % adjust to match QARTOD convention
    %IL(IL == 0) = 4; % this statement does NOT work because IL is a logical variable
    
    I = ones(size(IL,1),1);
    I(IL == 0) = 4;
    TUVqc.QC305 = I;  % there is no QARTOD location test but will call this QC305 for now
    clear I, IL;
    
    Uerr_test_parameters = {conf.QC.GDOP_OI.QC306_OI_Uerr.ErrorTypeStr conf.QC.GDOP_OI.QC306_OI_Uerr.ErrorVariableStr conf.QC.GDOP_OI.QC306_OI_Uerr.MaxError}';
    Verr_test_parameters = {conf.QC.GDOP_OI.QC307_OI_Verr.ErrorTypeStr conf.QC.GDOP_OI.QC307_OI_Verr.ErrorVariableStr conf.QC.GDOP_OI.QC307_OI_Verr.MaxError}';

    [~,I] = cleanTotals(TUVclean,conf.Totals.MaxTotSpeed,Uerr_test_parameters);  
    % I = Index of data points that were removed.  This will have a zero for
    %     good data, a 1 for data whose speed is above maxspd, a 2 for data
    %     violated first error condition, a 4 for data violating second error
    %     condition, etc.  For data violating multiple conditions, the result
    %     will be a sum of all the violated conditions.
    % b/c max speed already called index variable will not be a sum since both tests cannot fail
    I(I == 2) = 4;
    I(I == 0) = 1;
    TUVqc.QC306 = I;  % there is no QARTOD OI Uerr uncertainty test but will call this QC306 for now
    clear I;
    
    [~,I] = cleanTotals(TUVclean,conf.Totals.MaxTotSpeed,Verr_test_parameters);  
    % b/c max speed already called index variable will not be a sum since both tests cannot fail
    I(I == 2) = 4;
    I(I == 0) = 1;
    TUVqc.QC307 = I; % there is no QARTOD OI Uerr uncertainty test but will call this QC307 for now
    clear I;
    
    TUVqc.qc_operator_flag = ones(length(TUV.U),1).*2;
    
    TUVqc.PRIM = max([TUVqc.QC303,TUVqc.QC305,TUVqc.QC306,TUVqc.QC307,TUVqc.qc_operator_flag ],[],2);
    minqc = min([TUVqc.QC303,TUVqc.QC305,TUVqc.QC306,TUVqc.QC307,TUVqc.qc_operator_flag ],[],2);
    minqc1 = find(minqc == 1);
    prim2 = find(TUVqc.PRIM == 2);
    setto1 = intersect(minqc1, prim2); 
    if ~isempty(setto1)
      TUVqc.PRIM(setto1) = 1;  % if at least one other test has passed, set the primary flag to pass instead of not evaluated
    end
  
    TUV = TUVqc;
    TUV.MARACOOS_TUV_struct_version = 'undefined';
    
    TUVmetadata.conf = conf;
    TUVmetadata.radial_metadata = [];
    for rr = 1:size(RTUV,1)
      clear tmp;
      tmp = char(RTUV(rr).FileName);
      TUVmetadata.radial_metadata{rr} = tmp(max(strfind(tmp,'/'))+1:end);
    end
    TUVmetadata.radial_num_sites = size(RTUV,1);
    TUVmetadata.missingRadials = missingRadials;
    TUVmetadata.attributes = [];
    TUVmetadata.header = ['QCPrimaryFlagDefinition: Highest flag value of these test flags: ',conf.QC.OITotalPrimaryFlag_QC_Tests ,' (Note: Primary flag will be set to 2 only if ALL tests were not evaluated)'];
    

    % Save results
    [tdn, tfn] = datenum_to_directory_filename([conf.Totals.OI.BaseDir,...
        'mat/', conf.SystemType,'/', lower(pattern_field)],...
        dtime,...
        conf.Totals.OI.FilePrefix,...
        conf.Totals.FileSuffix,...
        conf.Totals.MonthFlag);
    
    tdn = tdn{1};
   
    if ~exist(tdn, 'dir')
      mkdir(tdn);
    end

    save(fullfile(tdn, tfn{1}),'TUV','RTUV', 'TUVmetadata' );
    
    
     [tdn_nc, ~] = datenum_to_directory_filename([conf.Totals.OI.BaseDir,...
        'nc/', conf.SystemType,'/realtime', lower(pattern_field)],...
        dtime,...
        conf.Totals.OI.FilePrefix,...
        '.nc',...
        0);
    tdn_nc = tdn_nc{1};
    if ~exist(tdn_nc, 'dir')
        mkdir(tdn_nc);
    end
 
   eval(['!/home/codaradm/operational_scripts/hfradarpy_scripts/bash/realtime_qc_totals_to_netcdf.sh ', fullfile(tdn, tfn{1}), ' ', tdn_nc ])



end
fprintf(1, '--------------------------------------------------------- \n');


