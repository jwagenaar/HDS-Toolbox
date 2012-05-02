function varargout = findHDS(obj, returnObj, varargin)
    % FINDHDS  Find objects in the database
    %
    %   RESULT = FINDHDS('ReturnClass', search syntax, ...)
    %   
    %   RESULT = FINDHDS('ReturnClass', OBJECT, search syntax,...)
    %
    % 
    % String variables will automatically be truncated to a maximum of 15
    % characters during the search. 

    global HDSsearchTables
    global HDSManagedData
    
    try

        % Check which index of HDSMManagedData should be searched.
        if isempty(HDSManagedData)
            throwAsCaller(MException('HDS:findhds','No database loaded in memory, unable to conduct search.'));
        elseif length(HDSManagedData) > 1
            if isa(varargin{1},'HDS')
                id = varargin{1}.objIds(1);
                treeIx = 0;
                for i = 1:length(HDSManagedData)
                    chk = any(HDSManagedData(i).objIds(1,:) == id);
                    if chk
                        treeIx = i;
                        break
                    end

                end
                if treeIx == 0
                    throwAsCaller(MException('HDS:findhds','FINDHDS: Unable to find database associated with supplied object.'))
                end
            else
                throwAsCaller(MException('HDS:findhds','FINDHDS: There are multiple paralell databases in memory. Include an object from the database you want to search through as the 3rd argument to the FINDHDS method to indicate the desired database.'));
            end
        else
            treeIx = 1;
        end

        if isempty(HDSsearchTables) 
            loadsearchtables(HDSManagedData(treeIx).basePath, HDSManagedData(treeIx).treeConst(1));
        elseif HDSManagedData(treeIx).treeConst(1) ~= HDSsearchTables.hostId
            loadsearchtables(HDSManagedData(treeIx).basePath, HDSManagedData(treeIx).treeConst(1));
        end

        % Form a struct based on the inputs describing the search.
        queryStruct = formQueryStruct(obj, returnObj, varargin);

        % If the returned object is equal to start object. Then there is only
        % one path, otherwise, find all possible paths
        
        
        if ~strcmp(obj, returnObj)
            startBooleans = false(length(queryStruct.metaClasses),1);
            
            % TODO: So this should work except that the isInverted should be a
            % vector for each path and not a single boolean.
            pathMatrix = zeros(200,200);
            isInverted = [];
            pathIndex = 1;
%             while ~all(startBooleans)
%                 startObj = queryStruct.metaClasses{find(~startBooleans,1)};
%                 [TpathMatrix pathOk TisInverted] = getPathMatrix(startObj, returnObj, 1);
%                 if ~pathOk
%                     throwAsCaller(MException('HDS:findhds','Error getting path, possible loop detected.'));
%                 else
%                     
%                 end
%                 
%                 startBooleans(find(~startBooleans,1)) = true;
%             end
            startObj = class(obj);
            [pathMatrix pathOk isInverted] = getPathMatrix(startObj, returnObj, 1);
            metaMatrix = findMetaMatrix(queryStruct, pathMatrix);
            possiblePaths = getPossiblePaths(metaMatrix, pathMatrix);
        else
            pathMatrix = find(strcmp(obj, {HDSsearchTables.mDataTree.classes.name}),1);
            metaMatrix = findMetaMatrix(queryStruct, pathMatrix);
            possiblePaths = 1;
            isInverted = false;
        end

        % If the path is inverted. Use trackback in end and swap the begin and
        % end object.
        trackBack = isInverted;
        if trackBack
            trackBackObj = returnObj;
            returnObj = class(obj);
        end

        % Iterate over all possible paths and find results. 
        returnId = find(strcmp(returnObj, {HDSsearchTables.mDataTree.classes.name}), 1);
        result = zeros(100,1);
        jx = 0;
        allResults = cell(length(possiblePaths), size(pathMatrix,2));
        for ii = 1: length(possiblePaths)
                        
            ix = 1;
            out = findMetaInTable(HDSsearchTables, pathMatrix(possiblePaths(ii),ix), [], metaMatrix, pathMatrix(possiblePaths(ii),:), 0);
            allResults{ii,ix} = out;
            if length(pathMatrix(1,:)) > 1
                while 1
                    ix = ix+1;
                    out = findMetaInTable(HDSsearchTables, pathMatrix(possiblePaths(ii),ix), out, metaMatrix, pathMatrix(possiblePaths(ii),:), 1);
                    allResults{ii,ix} = out;
                    if out.classId == returnId
                        break
                    end
                end            
            end

            if~isempty(out.classIdx)
                ll = length(out.classIdx);
                if jx + ll >= length(result)
                    result = [result ; zeros(max([100 ll]),1)]; %#ok<AGROW>
                end

                result(jx+1:jx+ll) = out.classIdx;
                jx = jx+ll;
            end

        end

        result = sort(result(1:jx));

        if trackBack
            % Track the requested class back from the returned class.
            trackBackId = find(strcmp(trackBackObj, {HDSsearchTables.mDataTree.classes.name}),1);

            % Create a table with all object Ids and parent IDs of all objects
            % on path except for last class. The pT table has some redundant
            % information as column 3 and 4 together are represented by column 2.
            objs = [allResults{pathMatrix(possiblePaths,:) ~= returnId & pathMatrix(possiblePaths,:)  ~=0}];
            resultObjs = [out.pId];

            cId = [objs.cId];
            classId = round((cId - floor(cId))*1e4);

            pT = [ [objs.pId]; cId ; classId; [objs.classIdx]]';
            [~, uPT, ~] = unique(pT(:,2));
            pT = pT(uPT,:);
            % ---

            % Now walk through table getting the parent ids until the parent
            % with the correct class is found. 
            tBObj = zeros(size(pT,1),1);
            ix = 0;
            while ~isempty(resultObjs)

                loc = find(ismembc(pT(:,2), sort(resultObjs)));

                correctClass = pT(loc, 3) == trackBackId;
                if any(correctClass)
                    ix2 = ix + sum(correctClass);
                    tBObj(ix + 1:ix2) = pT(loc(correctClass), 4);
                    ix = ix2;
                end
                resultObjs = pT(loc,1);  
            end

            % return the indeces of the correct class.
            result = tBObj(1:ix);        
        end

        % If the path is inverted. Use trackback in end and swap the begin and end object.
        if trackBack
            returnObj = trackBackObj;
            returnId  = trackBackId ;
        end

        % Now get the resulting objects from the query
        resultVector = cell(length(result),1);
        oids = HDSsearchTables.(returnObj).oid(result);
        pids = HDSsearchTables.(returnObj).pid(result);
        flinks = HDSsearchTables.(returnObj).fLink(result);
        for i = 1 : length(resultVector)
            aux = HDSsearchTables.mDataTree.fLinks(1:find(HDSsearchTables.mDataTree.fLinks(:,flinks(i)),1,'last'),flinks(i));
            resultVector{i} = [aux' pids(i) returnId oids(i)];
        end

        if ~isempty(resultVector)
            objs = HDS.getObjFromLoc(resultVector, treeIx);
        else
            objs = [];
        end

        if nargout <= 1
            varargout{1} = objs;
        elseif nargout == 2
            varargout{1} = objs;
            varargout{2} = struct('paths',pathMatrix(possiblePaths,:));
        else
            throw(MException('findHDS:toomanyoutputs','Too many output variables'));
        end
    catch ME
        keyboard
        throwAsCaller(ME)
    end

end