function varargout = hdsfind(returnClass, varargin)
    % HDSFIND  Searches the database for specific objects.
    %
    %   RESULT = HDSFIND('ReturnClass', search syntax, ...) Returns objects
    %       of class 'ReturnClass' that satisfy the query defined in the search
    %       syntax. It will use the default search options '-fba' (see below)
    %       which means that it returns objects if all conditions are met.
    %
    %       Example 1: RESULT = HDSFIND('Trial', 'date == 03/17/2011', ...
    %                            'subjectId == subject1'); 
    %
    %       will return all 'Trial' objects where:
    %       1)  The Trial object or the objects on a branch descending
    %           (children) from the Trial object match the query or,
    %       2)  The Trial object or the objects on a branch ascending
    %           (parents) from the Trial object match the query.
    %   
    %       Note that the method will not find objects where some of the query
    %       is match in ascending objects and other parts of the query is
    %       matched in descending objects. Use the ISMEMBER or UNIQUE
    %       method to combine results of two searches. 
    %
    %       Example 2: Consider the data-tree: Subject - Experiment - Trial
    %            
    %           AUX1 = HDSFIND('Experiment', 'Trial:nr == 4');   
    %           AUX2 = HDSFIND('Experiment', 'Subject:id == 1'); 
    %           RESULT1 = RESULT1(ISMEMBER(AUX1, AUX2));
    %           RESULT2 = UNIQUE([AUX1 AUX2]);
    %
    %       RESULT1 will return all Experiment objects that satisfy both
    %       search queries. RESULT2 will return all Experiment objects that
    %       satisfy either query. Note that the first query uses forward
    %       searching and the second query uses backwards searching.
    %   
    %   RESULT = HDSFIND('ReturnClass', OBJECT, search syntax,...) is used
    %       when multiple parallel databases are loaded into the MATLAB
    %       workspace. OBJECT should be any object from the database you want
    %       to search. It is only used to determine the specific database and
    %       does not further influence the search.
    %
    %   RESULT = HDSFIND('ReturnClass, '-searchOptions', ...) 
    %       uses custom options for the search. If omitted, the default options
    %       are used. Possible search-options are:
    %
    %       - f     : Forward searching of the database
    %       - b     : Backward searching of the database
    %       - a     : AND all search criteria
    %       - o     : OR all search criteria
    %
    %       The default searchOption is '-fba' which means that the database
    %       will search forwards and backwards and that all search criteria
    %       need to be met. Alternatives are '-fbo', '-fa', '-ba', '-fo', and
    %       '-bo'. 
    %
    %   RESULT = HDSFIND('ReturnClass', '-f', search syntax, ...) only searches
    %       paths where the returned objects are the last decendent of the
    %       searched paths. 
    %
    %   RESULT = HDSFIND('ReturnClass', '-b', search syntax, ...) only searches
    %       paths where the returned objects are the parent of the
    %       searched paths.
    %
    %
    %   == == == How does the search work == == ==
    %
    %   In order to understand the search function, you should realize that
    %   all data in the database is stored in a hierarchical data-tree
    %   structure. The search function searches for matches along the
    %   branches of the tree and returns objects of the requested class if
    %   all criteria are met on somewhere on the branch which includes the
    %   returned object. 
    %
    %   This means that if you have a data-tree that consists of three
    %   hierarchical classes: subject - experiment - trial, you can return
    %   all trial objects for a particular subject:
    %       For example: RESULT = HDSFIND('Trial','subjectId == subject1');
    %
    %   or you can return all subjects that have more than 10 trials.
    %       For example: RESULT = HDSFIND('Subject','trialNr > 10');
    %
    %   The latter can be interpreted as follows; If for a given Subject
    %   object, there is a Trial object where the property trialNr has a
    %   value > 10, then the search conditions are met and the Subject
    %   object is returned as a match.
    %
    %   Using the search-options variables in the syntax ('-fba') you can
    %   specify whether the search algorithm should include 'backward'
    %   results such as the previous example. 
    %
    %
    %   == == == Search Syntax == == ==
    %
    %   Syntax for the query can be specified in two different ways:
    %       1)  Defining a search term in string argument as in:   
    %           HDSFIND('class', 'property <comparator> argument') where the
    %           property, comparator and argument are all included in a single
    %           string. Multiple strings indicate multiple search criteria.
    %           For example:
    %               OUT = HDSFIND('className', 'channel == 5'); 
    %               OUT = HDSFIND('className', 'experimentName == Exp1');
    %
    %       2) Defining a search criterium using two arguments as in:
    %           OUT = HDSFIND('class', 'property <comparator>', argument) where
    %           the property and comparator are specified in a single string
    %           and the argument is defined separately.
    %           For example:
    %               OUT = HDSFIND('className', 'channel ==', 5);
    %               OUT = HDSFIND('className', 'channel <>',[1 2 3 4 5]);
    %               OUT = HDSFIND('className', 'experimentName ==', 'exp1');
    %
    %   The second syntax option is usefull in case the search criterium
    %   argument contains multiple values or is dynamically included.
    %
    %   When the argument can be represented by a number when using syntax
    %   option 1, the number will be automatically converted to a double.
    %
    %   More specific, you can also include which class should be checked
    %   for a property value. This can be usefull if multiple classes have
    %   the same property names. This is indicated in the following list:
    %
    %   - <'property comparator argument'> syntax:
    %       Use this syntax when you don't care which class contains the
    %       property that you match and if the argument is a single number
    %       or string.
    %           example: RESULT = HDSFIND('Trial', 'id == subject1');
    %
    %   - <class:property comparator argument> syntax:
    %       Use this syntax if you want to specify directly which class you
    %       want to check for the property and if the argument is a single
    %       number or string.
    %           example: RESULT = HDSFIND('Trial', 'Subject:id == subject1');
    %       
    %   - <class:property comparator> <argument> syntax:
    %       Use this syntax if you don't care which class contains the
    %       property that you match and if the argument is either assigned
    %       dynamically or contains multiple values.
    %           example: RESULT = HDSFIND('Trial', 'id <>',[ 1 2 3 4]);
    %
    %   - <class:property/property comparator> <argument> syntax:
    %       Use this syntax if you want to specify directly which class you
    %       want to check for the property and if the argument is either
    %       dynamically assigned or contains multiple variables.
    %           example: RESULT = HDSFIND('Trial', 'Subject:id <>', [ 1 2 3 4]);     
    %
    %
    %   == == == Comparator Symbols == == ==
    %
    %   Various comparator symbols can be specified to specify search
    %   parameters. The following symbols are allowed:
    %   
    %   '<'     : smaller than
    %   '<='    : smaller or equal to
    %   '=='    : equal to
    %   '~='    : not equal to
    %   '>='    : bigger or equal to
    %   '>'     : bigger than
    %   '<>'    : Include any (true if property contains any of arguments)
    %   '[]'    : Include all (true if property contains all of arguments)
    %
    %   The '<>' and '[]' symbols can be used on numeric and string
    %   arguments. The '[]' symbol results in the same behavior as the '=='
    %   symbol when used on strings.
    %
    %
    %   == == == Searching for Linked Objects == == ==
    %
    %   You can search match for linked objects that are added to an object
    %   using the ADDLINK method. To do so, use 'link' as the property name
    %   and include the linked object(s) as the argument. 
    %
    %   EXAMPLE: If you have two tree-branches: 
    %       subject - Experiment - Trial - RecordingChannel, and
    %       subject - ElectrodeArray - ElectrodeSite,
    %
    %       and you link each recordingChannel obejct to an ElectrodeSite
    %       object and an electrodeArray object using the ADDLINK method. 
    %
    %       In this scenario, you can search for all recordingChannel
    %       objects that were recorded on a specific ElectrodeSite using: 
    %       RESULT = HDSFIND('RecordingChannel','link ==', <ElectrodeSite object>);
    %
    %       or, all the trials that used a specific array:
    %       RESULT = HDSFIND('Trial','link ==', <ElectrodeArray object>);
    %       
    %
    %   See also: ISMEMBER, UNIQUE
    
    global HDSsearchTables
    
    try
        
        [searchOptions, treeIx, queryStruct] = hdsfindInputs(returnClass, varargin{:});
        
        % Find the Possible search paths and associated metaMatrix
        allClassNames  = {HDSsearchTables.mDataTree.classes.name};
        returnClassIdx = find(strcmp(returnClass, allClassNames ),1);
        returnClassIdxOriginal = returnClassIdx; % Save real returnclass index because returnclassidx can change
        returnClassOriginal = returnClass;
        
        % StartBooleans is used to check whether all paths from associated
        % class are included.
        startBooleans = true(length(queryStruct.metaClasses),1);
        startIndeces  = queryStruct.metaClassIdx;
        
        %PathMatrix is used to contain all possible paths
        pathMatrix = zeros(100,100);
        AvRowIdx   = 1;
        
        % Iterate over all possible paths. 
        while any(startBooleans)
            startClass = queryStruct.metaClasses{find(startBooleans,1)};
            pathMatrix_add = getPathMatrix(startClass, returnClass, searchOptions);
            
            pathMatrix(AvRowIdx:AvRowIdx+size(pathMatrix_add,1)-1, 1:size(pathMatrix_add,2)) = pathMatrix_add;
            AvRowIdx = AvRowIdx + size(pathMatrix_add,1);
            
            uniqueClassIdx = unique(pathMatrix_add);
            startBooleans = startBooleans & ~ismembc(startIndeces, uniqueClassIdx);
        end
        
        % metamatrix does something important
        metaMatrix = findMetaMatrix(queryStruct, pathMatrix);
        
        % possiblepaths return the indeces of possible paths --> rows of
        % pathmatrix.
        possiblePaths = getPossiblePaths(metaMatrix, pathMatrix);        
        
        % Iterate over paths and find resulting objects.

        allResults = cell(length(possiblePaths), size(pathMatrix,2));
        jx = 0;
        result = zeros(1000,2); % [index classIdx]

        for iPath = 1: length(possiblePaths)
            
            % Rearrange start and return objects for trackback.
            if pathMatrix(possiblePaths(iPath),1) == returnClassIdxOriginal && pathMatrix(possiblePaths(iPath),2) ~=0
                returnClass  = allClassNames{pathMatrix(possiblePaths(iPath),find(pathMatrix(possiblePaths(iPath),:),1,'last'))};  
                returnClassIdx = find(strcmp(returnClass, allClassNames ),1);
                trackBack  = true;
            else
                trackBack = false;
            end
            
            ix = 1;
            out = findMetaInTable(HDSsearchTables, pathMatrix(possiblePaths(iPath),ix), [], metaMatrix, pathMatrix(possiblePaths(iPath),:), 0);
            allResults{iPath,ix} = out;
            if pathMatrix(iPath,1) ~= returnClassIdx
                while 1
                    ix = ix+1;
                    out = findMetaInTable(HDSsearchTables, pathMatrix(possiblePaths(iPath),ix), out, metaMatrix, pathMatrix(possiblePaths(iPath),:), 1);
                    allResults{iPath,ix} = out;
                    if out.classId == returnClassIdx
                        break
                    end
                end            
            end
            
            % If trackBack than get back to original class Idx
            if trackBack
                resultObjs = [out.pId];
                curResults = allResults(possiblePaths(iPath),:);
                
                
                objs = [curResults{ pathMatrix(possiblePaths(iPath),:) ~= returnClassIdx & pathMatrix(possiblePaths(iPath),:)  ~=0}];

                cId = [objs.cId];
                classId = round((cId - floor(cId))*1e4);

                pT = [ [objs.pId]; cId ; classId; [objs.classIdx]]';
                [~, uPT, ~] = unique(pT(:,2));
                pT = pT(uPT,:);
                
                % Now walk through table getting the parent ids until the parent
                % with the correct class is found. 
                tBObj = zeros(size(pT,1),1);
                ix = 0;
                while ~isempty(resultObjs)

                    loc = find(ismembc(pT(:,2), sort(resultObjs)));

                    correctClass = pT(loc, 3) == returnClassIdxOriginal;
                    if any(correctClass)
                        ix2 = ix + sum(correctClass);
                        tBObj(ix + 1:ix2) = pT(loc(correctClass), 4);
                        ix = ix2;
                    end
                    resultObjs = pT(loc,1);  
                end
                
                % Create new out object
                out = struct('classId', returnClassIdxOriginal, 'classIdx' , tBObj(1:ix));
            end
            
            
            if~isempty(out.classIdx)
                ll = length(out.classIdx);
                if jx + ll >= length(result)
                    result = [result ; zeros(max([100 ll]),2)]; %#ok<AGROW>
                end

                result(jx+1:jx+ll,1) = out.classIdx;
                result(jx+1:jx+ll,2) = out.classId;
                jx = jx+ll;
            end
            
            
        end
        
        % truncate result and find unique indeces
        result = result(1:jx ,:);
        [~ , loc] = unique(result(:,1));
        result = result(loc,:);
        

        % Now get the resulting objects from the query
        resultVector = cell(length(result),1);
        oids = HDSsearchTables.(returnClassOriginal).oid(result);
        pids = HDSsearchTables.(returnClassOriginal).pid(result);
        flinks = HDSsearchTables.(returnClassOriginal).fLink(result);
        for i = 1 : length(resultVector)
            aux = HDSsearchTables.mDataTree.fLinks(1:find(HDSsearchTables.mDataTree.fLinks(:,flinks(i)),1,'last'),flinks(i));
            resultVector{i} = [aux' pids(i) returnClassIdxOriginal oids(i)];
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
        throwAsCaller(ME);
    end
end

function [searchOptions, treeIx, queryStruct] = hdsfindInputs(returnClass, varargin)
    % HDSFININPUTS  Analyzes the inputs of the HDSFIND method.
    
    global HDSManagedData
    
    % Find if there is an object specified indicating the database that
    % should be used.
    
    if length(varargin) > 1
    
        if isa(varargin{1}, 'HDS') || (isa(varargin{2}, 'HDS') && strcmp(varargin{1}, 'link'))
            if isa(varargin{1}, 'HDS'); objIndex = 1; else objIndex = 2; end;
            treeIx = varargin{objIndex}.treeNr;
            if ~treeIx
                [~, treeIx] = registerObjs(varargin{objIndex}, [],[],[], false);
            end
            varargin(objIndex) = [];

        elseif length(HDSManagedData) > 1
            throwAsCaller(MException('HDS:HDSFIND','HDSFIND: There are multiple databases loaded in the Matlab workspace; include an object from the specific database in the search syntax to specify which database to search (see help).'));
        else
            treeIx = 1;
        end

        % Find if there are specific search options that apply
        if varargin{1}(1) == '-' || (isa(varargin{1}, 'HDS') && varargin{2}(1) == '-')
            if varargin{1}(1) =='-'; optionIndex = 1; else optionIndex = 2; end;

            %Check if options make sense
            if sum(ismember('-fbao',varargin{optionIndex})) ~= length(varargin{optionIndex})
                throwAsCaller(MException('HDS:HDSFIND','HDSFIND: Incorrect definition of the search options.'));
            end

            searchOptions = ismember('fbo',varargin{optionIndex});
            varargin(optionIndex) = [];
        else
            searchOptions = [true true false];
        end

        if ~any(searchOptions(1:2))
            searchOptions(1:2) = true;
        end
        
    else
        searchOptions = [true true false];
        if length(HDSManagedData) > 1
            throwAsCaller(MException('HDS:HDSFIND','HDSFIND: There are multiple databases loaded in the Matlab workspace; include an object from the specific database in the search syntax to specify which database to search (see help).'));
        else
            treeIx = 1;
        end
    end
        
    
    %  -- -- Find the search structure. -- --

    % Form a struct based on the inputs describing the search.
    queryStruct = formQueryStruct(returnClass, treeIx, varargin);

end