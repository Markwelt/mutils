function [B,I]  = rmoutliers(A,varargin)
%RMOUTLIERS  Remove outliers from data (copied from Matlab R2021b)
%
%   B = RMOUTLIERS(A) detects and removes outliers from data. A can be a
%   vector, matrix, table, or timetable. If A is a vector, RMOUTLIERS
%   removes the entries detected as outliers. If A is a matrix or a table,
%   RMOUTLIERS detects outliers for each column and then removes the rows
%   containing outliers.
%
%   B = RMOUTLIERS(A,DIM) reduces the size of A along the dimension DIM.
%   Use DIM = 1 to remove rows and DIM = 2 to remove columns.
%   RMOUTLIERS(A,DIM) first calls ISOUTLIER(A) to detect outliers.
%
%   B = RMOUTLIERS(A,..., METHOD) specifies the method used to determine
%   outliers. METHOD must be one of the following: 'median' (default),
%   'mean', 'quartiles', 'grubbs', or 'gesd'.
%
%   B = RMOUTLIERS(A,..., 'percentiles', [LP UP]) detects as outliers all
%   elements which are below the lower percentile LP and above the upper
%   percentile UP. LP and UP must be scalars between 0 and 100 with
%   LP <= UP.
%
%   B = RMOUTLIERS(A,..., MOVMETHOD, WL) uses a moving window method to
%   determine contextual outliers instead of global outliers. MOVMETHOD can
%   be 'movmedian' or 'movmean'.
%
%   B = RMOUTLIERS(A,..., 'MinNumOutliers',N) removes rows (columns) that
%   contain at least N outliers. N must be an integer. By default, N = 1.
%
%   B = RMOUTLIERS(A,..., 'ThresholdFactor', P) modifies the outlier
%   detection thresholds by a factor P.
%
%   B = RMOUTLIERS(A,..., 'SamplePoints',X) specifies the sample points X
%   representing the location of the data in A for the moving window
%   methods 'movmedian' and 'movmean'. If the first input A is a table, X
%   can also specify a table variable in A.
%
%   B = RMOUTLIERS(A,..., 'MaxNumOutliers', MAXN) specifies the maximum
%   number of outliers for the 'gesd' method only.
%
%   B = RMOUTLIERS(A,..., 'DataVariables',DV) removes rows according to
%   outliers in table variables DV. The default is all table variables in
%   A. DV must be a table variable name, a cell array of table variable
%   names, a vector of table variable indices, a logical vector, a function
%   handle that returns a logical scalar (such as @isnumeric), or a table 
%   vartype subscript.
%
%   [B,I] = RMOUTLIERS(A,...) also returns a logical column (row) vector I
%   indicating which rows (columns) of A were removed.
%
%   Examples:
%
%     % Remove outliers from a vector
%       a = [1 2 1000 3 4 5]
%       b = rmoutliers(a)
%
%     % Remove only the rows which contain at least 2 outliers
%       A = [[1 2 1000 3 4 5]', [1 2 1000 3 4 1000]']
%       [B,removedRows] = rmoutliers(A,'MinNumOutliers',2)
%
%   See also ISOUTLIER, FILLOUTLIERS, RMMISSING, ISMISSING, FILLMISSING

%   Copyright 2018-2021 The MathWorks, Inc.

[B,I] = rmMissingOutliers('rmoutliers',A,varargin{:});
end
 

function [B,I] = rmMissingOutliers(funName,A,varargin)
% rmMissingOutliers Helper function for rmmissing and rmoutliers
%
%   FOR INTERNAL USE ONLY -- This feature is intentionally undocumented.
%   Its behavior may change, or it may be removed in a future release.
%
%   B - A after removing rows or columns
%   I - Colum(row) logical vector indicating removed rows(columns)
%

%   Copyright 2015-2018 The MathWorks, Inc.

opts = parseInputs(funName,A,varargin{:});

if ~opts.AisTable
    if ~ismatrix(A)
        issueError(funName,'NDArrays');
    end
    I = applyFun(funName,A,opts);
    I = cumputeIndex(I,opts.byRows,opts.minNum);
else
    if ~all(varfun(@ismatrix,A,'OutputFormat','uniform'))
        issueError(funName,'NDArrays');
    end
    if opts.byRows
        if opts.dataVarsProvided
            I = applyFun(funName,A(:,opts.dataVars),opts);
        else
            % Don't index into a table if we don't have to
            I = applyFun(funName,A,opts);
        end
        I = cumputeIndex(I,opts.byRows,opts.minNum);
        if isa(A,'timetable') && isequal(funName,'rmmissing')
            % Also remove the rows which correspond to missing RowTimes
            I = I | ismissing(A.Properties.RowTimes);
        end
    else
        I = false(1,width(A));
        for vj = opts.dataVars
            Ivj = applyFun(funName,A(:,vj),opts);
            I(vj) = cumputeIndex(Ivj,opts.byRows,opts.minNum);
        end
    end
end
B = reduceSize(A,I,opts.byRows);
end

%--------------------------------------------------------------------------
function I = applyFun(funName,A,opts)
% Find missing data or outliers
if isequal(funName,'rmmissing')
    I = ismissing(A);
else
    I = isoutlier(A,opts.isoutlierArgs{:});
end
end

%--------------------------------------------------------------------------
function I = cumputeIndex(I,byRows,minNum)
% Colum(row) logical vector indicating removed rows(columns)
I = sum(I,1 + byRows) >= minNum;
end

%--------------------------------------------------------------------------
function B = reduceSize(A,I,byRows)
% Keep non-missing
if byRows
    B = A(~I,:);
else
    B = A(:,~I);
end
end

%--------------------------------------------------------------------------
function opts = parseInputs(funName,A,varargin)
% Parse RMMISSING/RMOUTLIERS inputs

% We let the ISMISSING/ISOUTLIER calls error out for A of invalid type
opts.AisTable = matlab.internal.datatypes.istabular(A);

% Defaults
opts.minNum = 1;
opts.byRows = true;
opts.dataVarsProvided = false;
if ~opts.AisTable
    if isrow(A) && ~isscalar(A)
        opts.byRows = false;
    end
    opts.dataVars = NaN; % not supported for arrays
else
    opts.dataVars = 1:width(A);
end
if isequal(funName,'rmoutliers')
    opts.isoutlierArgs = {}; % arguments needed for ISOUTLIER computation
end

% Parse the outlier method, the trailing DIM and N-V pairs (including
% inputs which need to be forwarded to ISOUTLIER).
errorForDataVars = true;
opts = rmMissingOutliersVarargin(funName,A,opts,...
    errorForDataVars,varargin{:});
end

%--------------------------------------------------------------------------

function opts = rmMissingOutliersVarargin(funName,A,opts,...
    errorForDataVars,varargin)
% rmMissingOutliersVarargin Helper to parse DIM and N-V pairs for RMMISSING
% and RMOUTLIERS and inputs which need to be forwarded on to ISOUTLIER.
%
%   FOR INTERNAL USE ONLY -- This feature is intentionally undocumented.
%   Its behavior may change, or it may be removed in a future release.
%

%   Copyright 2018 The MathWorks, Inc.

numargs = numel(varargin);
if numargs == 0
    return
end
doOutliers = isequal(funName,'rmoutliers');

if doOutliers
    % Check if an outlier detection method has been specified
    input2 = varargin{1};
    offsetMethod = 0;
    if (ischar(input2) || isstring(input2))
        ind = matlab.internal.math.checkInputName(input2, ...
            {'median' 'mean' 'movmedian' 'movmean' 'percentiles' 'quartiles' 'grubbs' ...
            'gesd' 'SamplePoints' 'DataVariables' 'ThresholdFactor' ...
            'MinNumOutliers'});
        if sum(ind) ~= 1
            error(message('MATLAB:rmoutliers:SecondInputString'));
        end
        offsetMethod = any(ind(1:8)) + any(ind(3:5));
    end
    if offsetMethod == 0
        opts.isoutlierArgs = {'median'}; % Use the default method
    else
        % Pass the method on to isoutlier, including the invalid case of
        % isoutlier(A,movmethod)
        opts.isoutlierArgs = varargin(1:min(numargs,offsetMethod));
    end
    startNV = offsetMethod + 1;
else
    startNV = 1;
end

% Parse the DIM and N-V pairs
if numargs >= startNV
    dimIn = varargin{startNV};
    [opts,startNV] = getDim(funName,opts,startNV,dimIn,numargs,doOutliers);
    
    indNV = startNV:numargs;
    if rem(numel(indNV),2) ~= 0
        issueError(funName,'NameValuePairs');
    end
    
    if doOutliers
        minNumName = 'MinNumOutliers';
    else
        minNumName = 'MinNumMissing';
    end
    
    extraArgs = [];
    for k = indNV(1:2:end)
        opt = varargin{k};
        if matlab.internal.math.checkInputName(opt,minNumName)
            if doOutliers && matlab.internal.math.checkInputName(opt,'m')
                % 'm' is ambiguous due to 'MinNum' and 'MaxNum'
                validatestring('m',{'MinNumOutliers', 'MaxNumOutliers'});
            end
            minNum = varargin{k+1};
            if (~isnumeric(minNum) && ~islogical(minNum)) || ...
                    ~isscalar(minNum) || ~isreal(minNum) || ...
                    fix(minNum) ~= minNum || ~(minNum >= 0)
                issueError(funName,minNumName);
            end
            opts.minNum = minNum;
        elseif matlab.internal.math.checkInputName(opt,'DataVariables')
            opts.dataVarsProvided = true;
            dataVars = varargin{k+1};
            if errorForDataVars
                if opts.AisTable
                    dataVars = matlab.internal.math.checkDataVariables(...
                        A,dataVars,funName);
                else
                    issueError(funName,'DataVariablesArray');
                end
                opts.dataVars = dataVars;
            else
                % Just collect the data variables for later validation
                opts.dataVars = {opts.dataVars{:} dataVars}; %#ok<CCAT>
            end
        else
            extraArgs = [extraArgs k k+1]; %#ok<AGROW>
        end
    end
    % MinNum and DataVariables apply to both rmmissing and rmoutliers,
    % while everything else is a N-V for rmoutliers and isoutlier
    if doOutliers
        opts.isoutlierArgs = [opts.isoutlierArgs varargin(extraArgs)];
    elseif ~isempty(extraArgs)
        error(message('MATLAB:rmmissing:NameValueNames'));
    end
end
end

%--------------------------------------------------------------------------
function [opts,startNV] = getDim(funName,opts,startNV,dim,nargs,doOutliers)
% If an optional DIM is specified, then it must be the second input for
% RMMISSING and a trailing input for RMOUTLIERS (same as in ISOUTLIER):
if (doOutliers || nargs > 1) && (ischar(dim) || isstring(dim))
    %   rmmissing(A,N1,V1,N2,V2,...)
    %   rmoutliers(A,N1,V1,N2,V2,...)
    %   rmoutliers(A,method,N1,V1,N2,V2,...)
    %   rmoutliers(A,movmethod,window,N1,V1,N2,V2,...)
    % startNV = startNV; % N-V pairs can be the first entry of varargin
else
    %   rmmissing(A,DIM,N1,V1,N2,V2,...)
    %   rmoutliers(A,DIM,N1,V1,N2,V2,...)
    %   rmoutliers(A,method,DIM,N1,V1,N2,V2,...)
    %   rmoutliers(A,movmethod,window,DIM,N1,V1,N2,V2,...)
    startNV = startNV+1; % N-V pairs can be the second entry of varargin
    if (isnumeric(dim) || islogical(dim)) && isreal(dim) && isscalar(dim)
        if dim == 1
            opts.byRows = true;
        elseif dim == 2
            opts.byRows = false;
        else
            issueError(funName,'DimensionInvalid');
        end
    else
        issueError(funName,'DimensionInvalid');
    end
end
end

%--------------------------------------------------------------------------
function issueError(funName,errorId)
% Issue error from the correct error message catalog
error(message(['MATLAB:', funName, ':', errorId]));
end
