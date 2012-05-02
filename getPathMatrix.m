function varargout = getPathMatrix(fromClass, toClass, searchOptions, varargin)
    % GETPATHMATRIX  Method for determining the search path within database
    %
    % This methods is internally used to determine the different paths in
    % the database between objects. This is used to maximize efficiency for
    % the search algorithm such that it does not have to search through all
    % possible classes for the meta data.

    global HDSsearchTables
    
    % Init Path search
    direction       = 1;
    showProgress    = false;    % for debug, see path construction
    checkTree       = true;     % Check tree for loops
    loopDetect      = false;    % Check for loops
    hierarchyDetect = true;     % 
    isInverted      = false;    % Set to true if path is inverted.

    error           = [];
    maxActiveRows   = 10000;
    maxColnr        = 1000;
    
    switch nargin
        case 5
            direction = varargin{1};            
        case 6
            direction = varargin{1};
            showProgress = varargin{2};
        case 7
            direction = varargin{1};
            showProgress = varargin{2};
            checkTree = varargin{3};    
    end
    
    nrClasses = size(HDSsearchTables.mDataTree.links,1);

    toClassIdx   = find(strcmp(toClass, {HDSsearchTables.mDataTree.classes.name}), 1);
    fromClassIdx = find(strcmp(fromClass, {HDSsearchTables.mDataTree.classes.name}), 1);
    if isempty(toClassIdx) || isempty(fromClassIdx)
        throw(MException('getPathMatrix:IncorrectQuery', 'Incorrect query.'));
    end
    
    % If from class is the same as to class then return fromClassIdx
    if strcmp(fromClass,toClass)
        varargout = {fromClassIdx true false};
        return
    end
        
    
    [row ~] = find(HDSsearchTables.mDataTree.links);
    if ~any(HDSsearchTables.mDataTree.links(fromClassIdx,:)) && length(unique(row)) == nrClasses-1;
        hierarchyDetect = false;
    end
    
    if checkTree || hierarchyDetect
        % Generate hierarchy
        evalClass = ones(nrClasses,1);

        hierarchy = nan(nrClasses,1);
        while any(evalClass)
            ii = find(evalClass,1);
            childIds  = logical(HDSsearchTables.mDataTree.links(:,ii));
            parentIds = logical(HDSsearchTables.mDataTree.links(ii,:));

            if any(parentIds)
                hParent = max(hierarchy(parentIds));
            else
                hParent = max([0 min(hierarchy)])-1;
            end

            if isnan(hParent)
                hParent = 0;
            end

            if any(childIds)
                hChild = min(hierarchy(childIds));
            else
            	hChild = max(hierarchy)+1;
            end
            if isnan(hChild)
                hChild = max([0 max(hierarchy)])+1;
            end

            if hParent < hChild
                hierarchy(ii) = (hParent + hChild)/2;
            else
                hierarchy(ii) = hChild/2;
                evalClass(hierarchy(parentIds) >= hierarchy(ii))= true;
            end

            evalClass(ii) = false;
        end
    end
    
    % If checkTree is true then we should check for loops in the tree
    % structure. This is only necessary once per search.
%     if checkTree
%         % Check for loops
%         for ii =1 : nrClasses
%             parents = find(HDSsearchTables.mDataTree.links(ii,:));
%             if any(parents)
%                 if any(hierarchy(parents) > hierarchy(ii))
%                     % Found possible loop
%                     loopDetect = true;
% 
%                     whichParent = parents(find(hierarchy(parents) > hierarchy(ii), 1));
%                     Msg = sprintf(['Found possible loop in data structure; Cannot perform search.\n'...
%                         'Class %s can be a parent and a child of class %s\n'],...
%                         upper(HDSsearchTables.mDataTree.classNames{whichParent}),upper(HDSsearchTables.mDataTree.classNames{ii}));
%                     fprintf(2,Msg);
% 
%                     error = MException('getPathMatrix:LoopDetect',Msg);
% 
%                     % Change search to find loop
%                     toClassIdx = ii;
%                     fromClassIdx = whichParent;
%                     isInverted = true;
%                     break;
%                 end
%             end
%         end
%     end
    
    % Check In and Output classes
    if ~isInverted
        inId = find(strcmp(fromClass, {HDSsearchTables.mDataTree.classes.name}),1);
        outId = find(strcmp(toClass, {HDSsearchTables.mDataTree.classes.name}),1);
        if hierarchy(inId) > hierarchy(outId)
            toClassIdx = inId;
            fromClassIdx = outId;
            isInverted = true;
        end
    end

    % Get Paths
    nrows  = 1;
    iColnr = 1;
    
    pathMatrix      = zeros(50,20);
    pathMatrix(:,1) = toClassIdx;
    pathMatrixSize  = size(pathMatrix);
    
    availableRows   = 1:50;
    avStart = 1;
    avEnd   = 1;
    
    while 1
        iColnr = iColnr+1;
        
        % iRowActive is rows with possible branches
        % Rows in pathMatrix -> each row is possible path.
        iRowActive = find(pathMatrix(1:nrows, iColnr-1) ~= fromClassIdx & pathMatrix(1:nrows, iColnr-1) ~= 0);
        
        %break if iRowActive gets too large
        if length(iRowActive) > maxActiveRows || iColnr > maxColnr
            throw(MException('getPathMatrix:MaxActiveRow',sprintf('Datatree exceeded maximum (= %i) potential paths.',maxActiveRows)));
        end
        
        % find parents of previous nodes ( row = iRowActive(row) ) iRowActive_child in
        % mDataTree.links -> each iRowActive_child is classID and each iRowActive_parent is classID of
        % parent of class in iRowActive_child.
        activeObjects = pathMatrix(iRowActive, iColnr-1); %ClassIDs of Active Rows.
        [iRowActive_child iRowActive_parent] = find(HDSsearchTables.mDataTree.links(activeObjects, :));
        
        if showProgress
            display(sprintf('Col: %i -- Paths: %i',iColnr, length(iRowActive)));
        end
        
        % Check size of pathMatrix & grow if necessary
        if pathMatrixSize(1) < length(iRowActive_child) + nrows
            addmatrix = zeros(max([length(iRowActive_child) 100]) , pathMatrixSize(2));
            addmatrix(:,1) = toClassIdx;
            pathMatrix = [pathMatrix ; addmatrix]; %#ok<AGROW>
            pathMatrixSize = size(pathMatrix);
        end
        
        if pathMatrixSize(2) == iColnr
            pathMatrix = [pathMatrix zeros(pathMatrixSize(1), 20)]; %#ok<AGROW>
            pathMatrixSize = size(pathMatrix);
        end
        
        % Loop through paths
        changedIdx = zeros(nrows,1);
        n = 1;
        for ii = 1: length(iRowActive)                              % go through possible branches
            parentIds = iRowActive_parent(iRowActive_child == ii);  % the parents of the current object

            if ~isempty(parentIds)
                pathMatrix(iRowActive(ii), iColnr) = parentIds(1);
                changedIdx(n) = iRowActive(ii);
                n = n+1;
            
            
                for jj = 2: length(parentIds)    
                    
                    
                    if avStart ~= avEnd
                        avStart = avStart + 1;
                        if avStart > length(availableRows)
                            avStart = 1;
                        end

                        pathMatrix(availableRows(avStart), :) = pathMatrix(iRowActive(ii),:);
                        pathMatrix(availableRows(avStart), iColnr) = parentIds(jj);

                        changedIdx(n) = availableRows(avStart);
                        n = n+1;

                    else
                        nrows = nrows + 1;
                        pathMatrix(nrows, :) = pathMatrix(iRowActive(ii),:);
                        pathMatrix(nrows, iColnr) = parentIds(jj);

                        changedIdx(n) = nrows;
                        n = n+1;
                    end
                    
                end
            end
        end
        
        changedIdx = changedIdx(1:n-1);
                
        %Check length availableRows
        if avStart > avEnd
            if length(changedIdx) > avStart-avEnd
                l = length(changedIdx);
                availableRows = [availableRows zeros(1,max([100 l]))]; %#ok<AGROW>
                availableRows(l+1 : l+avEnd) = availableRows(avStart:l);
                avStart = l+1;
            end
        else
            if length(changedIdx) > length(availableRows)-(avEnd-avStart)
                l = length(changedIdx);
                availableRows = [availableRows zeros(1,max([100 l]))]; %#ok<AGROW>
            end
        end
        
        if ~isempty(changedIdx)
            if loopDetect
                for ii = 1: length(changedIdx)
                    if any(pathMatrix(changedIdx(ii), iColnr) == pathMatrix(changedIdx(ii),1:iColnr-1))
                        % Detect loops
                        avEnd = avEnd + 1;
                        if avEnd > length(availableRows)
                            avEnd = 1;
                        end

                        availableRows(avEnd) = changedIdx(ii);
                        pathMatrix(changedIdx(ii), iColnr) = 0;
                    end
                end
            elseif hierarchyDetect
                for ii = 1: length(changedIdx)    
                    if hierarchy(pathMatrix(changedIdx(ii), iColnr)) < hierarchy(fromClassIdx)
                        % Detect dead paths
                        avEnd = avEnd + 1;
                        if avEnd > length(availableRows)
                            avEnd = 1;
                        end

                        availableRows(avEnd) = changedIdx(ii);
                        pathMatrix(changedIdx(ii), iColnr) = 0;
                    end
                end
            end
        end

        % Break while loop when all paths are complete.
        if all(pathMatrix(:,iColnr)==0) 
            break
        end
        
    end
    
    % Get only complete paths
    [row col]  = find(pathMatrix == fromClassIdx);
    pathMatrix = pathMatrix(row,1:max(col));
    
    % Add row when LoopDetect
    if loopDetect
        pathMatrix = [repmat(fromClassIdx,size(pathMatrix,1),1) pathMatrix];
    end
    
    % Flip results
    if direction            
        for ii = 1: size(pathMatrix,1)
            v = pathMatrix(ii,pathMatrix(ii,:)>0);
            v = [v(end:-1:1) zeros(1,size(pathMatrix,2)-length(v))];
            pathMatrix(ii,:)=v;
        end
    end
    
    % Generate output
    if ~isempty(error)
        
        if nargout > 1
            if loopDetect
                varargout{2} = error;

                varargout{1} = struct('classes', [], 'loops', pathMatrix);
                uniqueClasses = unique(pathMatrix);
                uniqueClasses(uniqueClasses==0) = [];
                classNames = HDSsearchTables.mDataTree.classNames(uniqueClasses);
                classNames = classNames(:);
                classIds = num2cell(uniqueClasses);
                classIds = classIds(:);
                
                varargout{1}.classes = struct('className', classNames,'classId',classIds);
            else
                varargout{1} = pathMatrix;
                varargout{2} = error;
            end
        else
            varargout{1} = error;
        end
    else
        if nargout == 1
            varargout{1} = pathMatrix;
        elseif nargout == 2
            varargout{1} = pathMatrix;
            varargout{2} = true;
        elseif nargout == 3
            varargout{1} = pathMatrix;
            varargout{2} = true;
            varargout{3} = isInverted;
        end
    end
end