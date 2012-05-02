function hdsupdate(obj, varargin)
    %HDSUPDATE  Updates search lookup tables for given data tree.
    %
    %   HDSUPDATE(OBJ) Updates all search tables for the data tree
    %   associated with OBJ. Only changes to the data are updated.
    %
    %   HDSUPDATE(OBJ, 'all') Update all meta data including files that
    %   have not changed. This can potentially take a long time as it will
    %   browse through the entire database.
    %
    %   The HDSUPDATE method will generate an error report on the errors
    %   encountered while updating the search tables. The error report can
    %   be found in the 'HDSsearch' folder which is located at the root of
    %   the data tree.
    %
    %   In case an error is encountered, the HDSUPDATE method resolves this
    %   in the following matter:
    %   
    %   1) Incorrect class: The values in the property do not match the
    %   class of the initial value of the property as specified in the
    %   class definition. The exception is resolved by: ('boolean': false),
    %   ('numeric': NaN), ('char': '').
    %
    %   2) Incorrect size: The values in the property do not have the
    %   correct dimensions. For 'booleans' and 'numeric' classes, the
    %   property should contain a single value. For 'char' the property
    %   should contain a 1xn char array. The exception is resolved by:
    %   ('boolean': first value in property), ('numeric': NaN), 
    %   ('char': '').
    %
    %   3) Empty index: The property is empty. No value is present in the
    %   property. The exception is resolved by: ('numeric': NaN), 
    %   ('char': '').

    % Copyright (c) 2012, J.B.Wagenaar
    % This source file is subject to version 3 of the GPL license, 
    % that is bundled with this package in the file LICENSE, and is 
    % available online at http://www.gnu.org/licenses/gpl.txt
    %
    % This source file can be linked to GPL-incompatible facilities, 
    % produced or made available by MathWorks, Inc.
            
    global HDSManagedData 
    global HDSsearchTables
    

    if length(obj)>1
        throwAsCaller(MException('HDS:hdsupdate','Method not defined for arrays of objects'));
    end
    
    % Clear HDSsearchtables because it should be updated.
    HDSsearchTables = [];
    
    % Register object if not previously registered.
    treeId = obj.treeNr;
    if ~treeId
        [~, treeId] = registerObjs(obj);
    end
    
    basePath = HDSManagedData(treeId).basePath;
    
    % Check additional input arguments.
    if nargin>1
        if any(strcmp(varargin{1},{'-all' 'all'}))
            includeAll = true;
        else
            throwAsCaller(MException('HDS:hdsupdatesearch','Incorrect input argument for HDSUPDATE method.'));
        end
    else
        includeAll = false;
    end
            
    % Load the HDSconfig file, this should be present when the data tree is
    % correctly saved.
    if exist(fullfile(basePath,'HDSconfig.mat'),'file')
        config = load(fullfile(basePath, 'HDSconfig.mat'));
    else
        throwAsCaller(MException('HDS:hdsupdate','Unable to find the ''HDSconfig.mat'' file.'));
    end
    
    %Find HDSsearch folder for lookup files or create folder
    searchPath = fullfile(basePath,'HDSsearch');
    if ~exist(searchPath,'dir')
        mkdir(basePath,'HDSsearch');
    elseif includeAll
        rmdir(searchPath,'s');
        mkdir(basePath,'HDSsearch');
    end
    
    % Find the current Tree structure and load the matrix with folder links
    % from disk if it exists.
    treeStruct = hdstreestruct(treeId);
    if exist(fullfile(searchPath, 'HDS_Treestruct.mat'),'file')
        allLinks = load(fullfile(searchPath, 'HDS_Treestruct.mat'),'fLinks');
        allLinks = allLinks.fLinks;
    else
        allLinks = zeros(6,1);
    end
    
    try 
        %Define a structure for the Update-report and create errorLog file.
        report = struct('nrFolders',0,'nrFiles',0,'nrAddedObjects',0);
        if includeAll
            errorfid = fopen(fullfile(searchPath,'errorlog.txt'),'w');
        else
            errorfid = fopen(fullfile(searchPath,'errorlog.txt'),'a');
        end

        % Write current date to errorlog.
        fprintf(errorfid, [datestr(clock) '\n']);
        
      % -- Iterate over classes of the data tree and update lookup table for each class.
        warning('off', 'MATLAB:class:mustReturnObject');
        classes = config.classes(1: (find(cellfun(@isempty,config.classes),1)-1) );   
        allClassIds = cell(length(classes), 2);
        linkedClassesInClass = cell(length(classes),1);
        
        for iClass = 1: length(classes)

            % Create testObject to find meta Properties. The create cell array
            % with the class of the metaProps; should be double, single, logical or str.
			hdspreventreg(true);
            tobj = eval(classes{iClass});
			hdspreventreg(false);
            metaProps = tobj.metaProps;
            metaPropClass = cell(length(metaProps),1);
            for iP = 1:length(metaProps)
                cl = class(tobj.(metaProps{iP}));
                if isempty(cl); cl = 'double';end;
                metaPropClass{iP} = cl;
            end

            % Check search file and create new file if it does not exist. Load
            % file it does exist and assign dataTable variable.
            searchFname = sprintf('S%s.mat',classes{iClass});
            if exist(fullfile(searchPath,searchFname),'file')
                dataTable = load(fullfile(searchPath, searchFname));
                optionstr = 'Updating';
            else
               % Define data table:
               % pid = parentId
               % oid = objectId
               % parentId = class Id of parent
               % parentIdx = index of parent in lookup table
               % comboId = parentId*0.0001 + parentIdx
               % fLink = []
               % linkId = []
               % linkClassId = []

               dataTable = struct('pid',[],'oid',[],'parentId',single([]),'parentIdx',[],'comboId',[],'fLink',[],'linkId',[],'linkClassId',[]);
               save(fullfile(searchPath, searchFname),'-struct','dataTable');

               % Create properties for metaData
               for iM = 1: length(metaProps)
                   switch metaPropClass{iM}
                       case 'double'
                           dataTable.(metaProps{iM}) = [];
                       case 'char'
                           dataTable.(metaProps{iM}) = '';
                       case 'logical'
                           dataTable.(metaProps{iM}) = false(0);
                       case 'single'
                           dataTable.(metaProps{iM}) = zeros(0, 'single');
                       case 'cell'
                           dataTable.([metaProps{iM} '_i']) = [];
                           dataTable.(metaProps{iM}) = '';
                   end
               end
               optionstr = 'Creating';
            end

            display(sprintf('(%i/%i) %s table for %s ', iClass, length(classes), optionstr ,upper(classes{iClass})));
            drawnow; % To force update the command window.

            % Get all the folders in which objects of the current class have
            % changed since the last update.
            changeFolders = getChangedFolders(basePath, classes{iClass}, includeAll, treeStruct);

            % Update report
            report.nrFolders = report.nrFolders + size(changeFolders,1);
            report.nrFiles   = report.nrFiles + sum(cellfun(@length, changeFolders(:,2)));

            % Sort the object IDS in the current dataTable
            [sortedIds, sortIdx] = sort(dataTable.oid);

            % Iterate over all changed folders,files and repopulate lookup table.
            for iFolder = 1: size(changeFolders,1)

                
                % Get FolderLink
                binFname = sprintf('%s.bin',classes{iClass} );
                fid = fopen(fullfile(changeFolders{iFolder,1}, binFname),'r+');
                data = fread(fid,'*uint32');
                Flink = data(4: 3+data(3));
                if isempty(Flink); Flink = single(0); end;
                fclose(fid);

                % Check if FolderLink exists in allLinks and get the index. If
                % not, add the folderLink and set the index.
                allLinksSize = size(allLinks);

                % increase number of rows if Flink length is too large for matrix.
                if length(Flink) > allLinksSize(1)
                    allLinks = [allLinks ; zeros(length(Flink)-allLinksSize(1) + 2, allLinksSize(2) )]; %#ok<AGROW> % add two extra rows + length FLink.
                end

                
                % Check if Flink already exists.
                allLinksSize = size(allLinks);
                pFlink = double([Flink ; zeros(allLinksSize(1)-length(Flink),1)]); 
                pFlinkArray = pFlink(:, ones(allLinksSize(2),1));
                FlinkIndex = find(~any(allLinks - double(pFlinkArray,1)),1);

                % If it does not exist add to array.
                if isempty(FlinkIndex)
                    % Add colums if needed
                    allLinksSize = size(allLinks);
                    firstAvIndex = find(allLinks(1,:),1,'last') +1;
                    if isempty(firstAvIndex); firstAvIndex=2;end
                    if firstAvIndex > allLinksSize(2)
                        allLinks = [allLinks  zeros(allLinksSize(1), 10 )]; %#ok<AGROW> % add ten extra colums.
                    end

                    allLinks(:, firstAvIndex) = pFlink;
                    FlinkIndex = firstAvIndex;                    
                end

                % Iterates over files in Folder
                for iFile = 1: length(changeFolders{iFolder,2})

                    % If we are looking at main object --> iclass ==1 change file.
                    if strcmp(changeFolders{iFolder,1}, basePath)
                        curFile = fullfile(changeFolders{iFolder,1},sprintf('%s.mat',classes{iClass}));
                        objs = load(curFile);
                    else
                        curFile = fullfile(changeFolders{iFolder,1},sprintf('%s_%i.mat',classes{iClass},changeFolders{iFolder,2}(iFile)));
                        objs = load(curFile);
                    end


                    dataNames = fieldnames(objs);
                    out = objs.(dataNames{1});
                    out = out(ones(length(dataNames),1));

                    % Create array of new objects.
                    for i = 2:length(dataNames)
                        out(i) = objs.(dataNames{i}); 
                    end

                    objIds = {out.objIds};
                    objIds = [objIds{:}];

                    % Check if the objects already exist in the lookup table
                    objLoc = ismembc2(objIds(1,:), sortedIds);
                    objLoc(logical(objLoc)) = sortIdx(objLoc(logical(objLoc))); %Unsort objects that were found.

                    nrAddedObj = sum(~objLoc);
                    report.nrAddedObjects = report.nrAddedObjects + nrAddedObj;

                    % Create new indeces for unfound objects.                
                    objLoc(~objLoc) = (length(dataTable.oid)+1) : (length(dataTable.oid) + nrAddedObj);

                    dataTable.oid = [dataTable.oid zeros(1, nrAddedObj)];
                    dataTable.pid = [dataTable.pid zeros(1, nrAddedObj)];
                    dataTable.fLink = [dataTable.fLink zeros(1, nrAddedObj,'single')];
                    
                    for iM = 1:length(metaProps)
                        switch metaPropClass{iM}
                            case 'double'
                                dataTable.(metaProps{iM}) = [dataTable.(metaProps{iM}) zeros(1, nrAddedObj )];
                            case 'char'
                                space = ' ';
                                dataTable.(metaProps{iM}) = [dataTable.(metaProps{iM}) space(ones(size(dataTable.(metaProps{iM}),1), nrAddedObj ))];
                            case 'logical'
                                dataTable.(metaProps{iM}) = [dataTable.(metaProps{iM}) false(1, nrAddedObj )];
                            case 'cell'
                                dataTable.([metaProps{iM} '_i']) = [dataTable.([metaProps{iM} '_i']) zeros(1, nrAddedObj )];
                        end
                    end

                    % Populate OID
                    dataTable.oid(objLoc) = objIds(1,:);

                    % Populate PID
                    dataTable.pid(objLoc) = objIds(3,:);
                    
                    % Populate fLink
                    dataTable.fLink(objLoc) = single(FlinkIndex);
                    
                    % Populate linkClass                    
                    linkClasses = zeros(length(out), 0);
                    for i=1:length(out)
                        aux = double(out(i).linkPropIds(1, out(i).linkPropIds(2,:)==uint32(0)));
                        if length(aux) > size(linkClasses,2)
                            linkClasses = [linkClasses zeros(length(out), length(aux) - size(linkClasses,2))]; %#ok<AGROW>
                        end
                        linkClasses(i,1:length(aux)) = aux;
                    end
                    if ~isempty(linkClasses)
                        if size(dataTable.linkId,1) < size(linkClasses,2)
                            dataTable.linkClassId = [dataTable.linkClassId ; zeros(size(linkClasses,2)-size(dataTable.linkClassId,1), size(dataTable.linkClassId,2))];
                        end

                        dataTable.linkClassId(1:size(linkClasses,2), objLoc) = linkClasses';
                    end
                    
                    % Populate linkIds
                    linkIds = zeros(length(out), 0);
                    for i=1:length(out)
                        aux = out(i).linkIds(1, out(i).linkIds(3,:)~=0);
                        if length(aux) > size(linkClasses,2)
                            linkIds = [linkIds zeros(length(out), length(aux) - size(linkIds,2))]; %#ok<AGROW>
                        end
                        linkIds(i,1:length(aux)) = aux;
                    end
                    if ~isempty(linkIds)
                        if size(dataTable.linkId,1) < size(linkIds,2)
                            dataTable.linkId = [dataTable.linkId ; zeros(size(linkIds,2)-size(dataTable.linkId,1), size(dataTable.linkId,2))];
                        end

                        dataTable.linkId(1:size(linkIds,2), objLoc) = linkIds';
                    end
                    
                    % Populate PARENTID --> classId of parent
                    if strcmp(changeFolders{iFolder,1}, basePath)
                        dataTable.parentId(objLoc) = single(0);
                    else
                        [~, parentClassStr] = fileparts(changeFolders{iFolder,1});
                        parentClassStr = regexp(parentClassStr,'_','split');
                        parentClassStr = parentClassStr{1};
                        dataTable.parentId(objLoc) = single(find(strcmp(parentClassStr,classes),1));
                    end

                    % populate PARENTIDX, COMBOID with zeros
                    dataTable.parentIdx = zeros(1, length(dataTable.parentId), 'single');
                    dataTable.comboId = zeros(1, length(dataTable.parentId), 'single');
                    
                    
                    
                    

                  % -- Now add the values for the meta Properties --

                    % Currently, we only allow for doubles and character array.
                    % Strings will be transformed in a character array and
                    % srtings in properties will be truncated at 15 characters.

                    indexString = ''; %Init indexString
                    for iM = 1:length(metaProps)

                        if length(out) == 1
                            metaValues = out.(metaProps{iM});
                        else
                            metaValues = {out.(metaProps{iM})};
                        end

                        % Check class of metaValues; Values that do not have the same class are omitted.
                        if length(out) == 1
                           errorClass = ~isa(metaValues, metaPropClass{iM});
                        else
                           errorClass = ~cellfun('isclass', metaValues, metaPropClass{iM}); 
                        end

                        switch metaPropClass{iM}
                            case {'double' 'single'}

                                if length(out)==1
                                    % replace empty indeces with NaN
                                    if isempty(metaValues);
                                        metaValues = NaN;
                                        errorEmpty = true;
                                    else
                                        errorEmpty = false;
                                    end;

                                    % check length == 1; Values with different length are omitted.
                                    vallength = length(metaValues);
                                    errorLength = vallength~=1;

                                    useVals = ~errorClass & ~errorLength;
                                    metaValues(~useVals) = NaN;
                                    dataTable.(metaProps{iM})(objLoc) = metaValues;
                                else
                                    % replace empty indeces with NaN
                                    errorEmpty = cellfun('isempty',metaValues);
                                    metaValues(errorEmpty) = {NaN};

                                    % check length == 1; Values with different length are omitted.
                                    vallength = cellfun('length', metaValues);
                                    errorLength = vallength~=1;

                                    useVals = ~errorClass & ~errorLength;
                                    metaValues(~useVals) = {NaN};
                                    dataTable.(metaProps{iM})(objLoc) = [metaValues{:}];
                                end                                
                            case 'char'             
                                % In case of char, we check the length of the longest string. If < 15 we
                                % padd the existing array with spaces and add srings to array. This is
                                % better than keeping it as a cell array because a cell array takes way
                                % more memory and is ridiculousluy slow to load.

                                if length(out) == 1
                                    errorEmpty = isempty(metaValues);
                                    errorLength = size(metaValues,1)~=1;
                                    useVals = ~errorClass & ~errorLength;
                                    if ~useVals 
                                        metaValues = ' ';
                                    end
                                else
                                    errorEmpty = cellfun('isempty',metaValues);
                                    errorLength = cellfun('size',metaValues,1)~=1;
                                    useVals = ~errorClass & ~errorLength;
                                    metaValues(~useVals) = {' '};  % Set all metaValues that are incorrect to an empty string.
                                end

                                Values = char(metaValues)';

                                % Truncate the Values at a maximum of 15 characters or pad with spaces if needed.
                                sT = size(dataTable.(metaProps{iM}));
                                if size(Values,1) > 15
                                    Values = Values(1:15,:);
                                elseif size(Values,1) < sT(1)
                                    space = ' ';
                                    Values = [Values ; space(ones(sT(1) - size(Values,1),size(Values,2)))]; %#ok<AGROW>
                                end

                                % Pad dataTable with spaces if new inputs are longer than current table.
                                curLength = size(dataTable.(metaProps{iM}),1);
                                if curLength < size(Values,1)
                                    space = ' ';
                                    dataTable.(metaProps{iM}) = [dataTable.(metaProps{iM}) ; space(ones(size(Values,1)-curLength,size(dataTable.(metaProps{iM}),2)))];
                                end

                                try
                                dataTable.(metaProps{iM})(:,objLoc) = Values;
                                catch ME
                                    ME
                                    keyboard
                                end
                            case 'logical'          
                                % Impossible to have empty booleans so don't have to check for that.
                                errorEmpty  = false(length(metaValues),1);

                                if length(out) == 1      
                                    % Set value to false if errorClass
                                    if errorClass
                                        metaValues = false;
                                    end

                                    % check length == 1; Values with different length will only use first index.
                                    vallength   = length(metaValues);
                                    errorLength = vallength~=1;

                                    if errorLength
                                        metaValues = metaValues(1);
                                    end

                                    dataTable.(metaProps{iM})(objLoc) = metaValues;
                                else
                                    % Set indeces with incorrect class to false.
                                    if any(errorClass)
                                        metaValues(errorClass) = {false};
                                    end

                                    % check length == 1; Values with different length will only use first index.
                                    vallength   = cellfun('length', metaValues);
                                    errorLength = vallength~=1;

                                    if any(errorLength)
                                        metaValues = cellfun(@(x) x(1), metaValues,'uniformOutput',false);
                                    end

                                    dataTable.(metaProps{iM})(objLoc) = [metaValues{:}];
                                end 
                            case 'cell'             

                                if length(out) == 1
                                    errorEmpty = isempty(metaValues);
                                    errorLength = ~isvector(metaValues);
                                    useVals = ~errorClass & ~errorLength;

                                    % Cast as string if value is not string. If not all strings, try to cast as
                                    % strings. This works for boolean, string and numerics. Not for cells.
                                    if useVals && ~iscellstr(metaValues)
                                        if ~any(cellfun(@(x) iscell(x),metaValues))
                                            metaValues = cellfun(@(x) num2str(x), metaValues, 'uniformOutput',false);
                                        else
                                            useVals = false;
                                            errorClass = true;
                                        end
                                    end

                                else
                                    errorEmpty = cellfun('isempty',metaValues);
                                    errorLength = ~cellfun(@(x) isvector(x), metaValues);
                                    useVals = ~errorClass & ~errorLength;



                                    if ~all(cellfun(@(x) iscellstr(x) || ischar(x), metaValues))
                                        for i =1: length(out)
                                            % Cast as string if value is not string. If not all strings, try to cast as
                                            % strings. This works for boolean, string and numerics. Not for cells.
                                            if useVals(i) && ~iscellstr(metaValues{i})
                                                if ~any(cellfun(@(x) iscell(x), metaValues{i}))
                                                    metaValues{i} = cellfun(@(x) num2str(x), metaValues{i}, 'uniformOutput',false);
                                                else
                                                    useVals(i) = false;
                                                    errorClass(i) = true;
                                                end
                                            end
                                        end
                                    end

                                end

                                % At this point all the useVals indeces of
                                % metaValues are either strings or cell arrays
                                % of strings.
                                keyboard

                                % Concatenate strings in metaValues, if error
                                % it is because vector is column instead of
                                % row. then rotate cell array.
                                try
                                    metaValuesCC = [metaValues{:}];
                                catch

                                end    
                        end

                      % -- Write error information to the error log.
                        if any(errorLength) || any(errorClass) || any(errorEmpty)

                            % Find the indexString of the current file. Only do
                            % this once per file so first check if we computed
                            % the string already.
                            if isempty(indexString)
                                curFolder = changeFolders{iFolder,1}(length(basePath)+2:end);
                                indexString = '';
                                nF = [0 find(curFolder == filesep) length(curFolder)+1];
                                for ii = 1:length(nF)-1
                                    aux = curFolder((nF(ii)+1):(nF(ii+1)-1));
                                    splitStr = find(aux == '_',1);
                                    if isempty(splitStr)
                                        indexString = [indexString aux(1:min([3 length(aux)]))]; %#ok<AGROW>
                                    else
                                        indexString = [indexString sprintf('(%s).%s',aux((splitStr+1):end), aux(1:min([3 length(aux(1:splitStr-1))])))]; %#ok<AGROW>
                                    end
                                end
                                indexString = [indexString sprintf('(%i).%s',changeFolders{iFolder,2}(iFile), classes{iClass}(1:min([3 length(aux)])))]; %#ok<AGROW>

                                % Get order of indeces
                                fieldn = fieldnames(objs);
                                fieldidx = cellfun(@(x) str2double(x(2:end)), fieldn);
                            end

                            aux = find(errorEmpty);
                            if ~isempty(aux)
                                indeces = sprintf('%i ',sort(fieldidx(aux)));
                                fprintf(errorfid, 'Empty indeces: %s(%s).%s\n',indexString, indeces(1:end-1), metaProps{iM}); 
                            end

                            aux = find(errorClass);
                            if ~isempty(aux)
                                indeces = sprintf('%i ',sort(fieldidx(aux)));
                                fprintf(errorfid, 'Incorrect class: %s(%s).%s\n',indexString, indeces(1:end-1), metaProps{iM}); 
                            end

                            aux = find(errorLength);
                            if ~isempty(aux)
                                indeces = sprintf('%i ',sort(fieldidx(aux)));
                                fprintf(errorfid, 'Incorrect size: %s(%s).%s\n',indexString, indeces(1:end-1), metaProps{iM}); 
                            end
                        end

                    end

                    clear out
                end

              % -- Now open .bin file in folder and set the changedBooleans to
              % false. 
                try
                    binFname = sprintf('%s.bin',classes{iClass} );
                    fid = fopen(fullfile(changeFolders{iFolder,1}, binFname),'r+');
                    data = fread(fid,'*uint32');

                    % Set changedBooleans to zero.
                    data((5+data(3)): (4 + data(3) + data(4+data(3)))) = uint32(0);
                    frewind(fid);
                    fwrite(fid, data,'uint32');
                    fclose(fid);



                catch ME
                    fclose(fid);
                    rethrow(ME);
                end

            end

            % Iterate over de string metaProps and remove trailing spaces
            % ecept for the first row.

            strProps = find(strcmp('char', metaPropClass));
            for iS = 1:length(strProps)
                s = size(dataTable.(metaProps{strProps(iS)}));
                space = ' ';
                spaceArray = space(ones(s(1),s(2)));
                emptyRows = sum(dataTable.(metaProps{strProps(iS)}) - spaceArray,2) == 0;
                emptyRows(1) = false; % Make first row stay.
                dataTable.(metaProps{strProps(iS)})(emptyRows,:) = [];
            end

            linkedClassesInClass{iClass} = unique(dataTable.linkClassId);
            sparseLink = sparse(false(length(linkedClassesInClass)));
            for i = 1: length(linkedClassesInClass)
                sparseLink(i,linkedClassesInClass{i}) = true;
            end
            
            [allClassIds{iClass,1}, allClassIds{iClass,2}] = sort(dataTable.oid);
            save(fullfile(searchPath, searchFname),'-struct','dataTable');
            clear dataTable

        end

      % -- Save the Tree Structure --
        treeStruct.fLinks = allLinks;
        treeStruct.linkClass = sparseLink;
        save(fullfile(searchPath, 'HDS_Treestruct.mat'),'-struct','treeStruct');

      % -- Now, reopen all search files and populate the Parentidx --
        for iClass = 1: length(classes)
            searchFname = sprintf('S%s.mat',classes{iClass});

            parentTable = load(fullfile(searchPath, searchFname),'parentIdx','pid','parentId','comboId');

            emptyIdx = ~parentTable.parentIdx;
            sPID     = parentTable.pid(emptyIdx);

            % address problem with main object
            parentTable.parentIdx(sPID==0) = single(0);
            emptyIdx(sPID==0) = [];
            sPID(sPID==0)=[];

            uniqueParentIds = unique(parentTable.parentId(emptyIdx));
            for iP = 1:length(uniqueParentIds)
                pidx = ismembc2(sPID, allClassIds{uniqueParentIds(iP),1});
                if any(pidx==0)
                    keyboard
                    throwAsCaller(MException('HDS:update','Unable to find parent objects in table. Please rerun the HDSUPDATE method with the ''all'' attribute to correct this error.'));
                end
                parentTable.parentIdx(emptyIdx) = allClassIds{uniqueParentIds(iP),2}(pidx);
            end

            parentTable.comboId(emptyIdx) = parentTable.parentId(emptyIdx)*1e-4 + parentTable.parentIdx(emptyIdx);
%             parentTable.comboId(emptyIdx) = parentTable.parentId(emptyIdx)*1e10 + parentTable.parentIdx(emptyIdx);

            % Save the parent indeces to disk.
            parentIdx = parentTable.parentIdx; %#ok<NASGU>
            comboId  = parentTable.comboId; %#ok<NASGU>
            save(fullfile(searchPath,searchFname),'parentIdx','comboId','-append');

            clear parentIdx
        end

        HDS.displaymessage('-- -- -- -- -- --',2,'\n','');
        fprintf('HDSUPDATE completed updating the search tables.\n');
        fprintf('Number of updated folders :    %i \n',report.nrFolders);
        fprintf('Number of updated files   :    %i \n',report.nrFiles);
        fprintf('Number of added objects   :    %i \n',report.nrAddedObjects);
        HDS.displaymessage('-- -- -- -- -- --',2,'','\n');

        % Close the errorlog file
        fclose(errorfid); 
    catch ME
        fclose(errorfid);
        try fclose(fid); catch;end; %#ok<CTCH>
        rethrow(ME);
    end
end

function folders = getChangedFolders(basePath, className, includeAll, treeStruct)
    %GETCHANGEDFOLDERS  returns a nx2 cell array where the first column
    %are the paths to folders in which objects of 'className' have changed.
    %The second column indicates which class files have changed.
    
    returnFolders = cell(100,2);
    retIx       = 1;
    
    folders     = cell(100,2);
    folders(1,:)  = {basePath true};
    ix          = 1; % active row.
    ix2         = 2; % first available row
    
    clsnms = {treeStruct.classes.name};
    
    classIndex  = strcmp(className, clsnms);
    
    posParents  = clsnms(treeStruct.links(classIndex,:));
    
    while 1
      % -- Check if class exists in current dir and if so, if it has changed.
        if folders{ix,2}
            classPath = fullfile(folders{ix,1},sprintf('%s.bin',className));

            if exist(classPath,'file')

                fid = fopen(classPath);
                aux = fread(fid, 3,'*uint32');
                
                aux2 = fread(fid, aux(3)+1,'*uint32');
                nrFiles = aux2(end);
                
                ClassFileChanged = find(logical(fread(fid, nrFiles,'*uint32')));
                filesExist       = find(logical(fread(fid, nrFiles,'*uint32')));
                fclose(fid);

                if includeAll
                    % *Include all folders that contain objects of the current class.
                    returnFolders{retIx, 1} = folders{ix,1};
                    returnFolders{retIx, 2} = filesExist;
                    retIx = retIx+1;
                else
                    % *Only include folder if objects have changed.

                    % If class changed in folder than add folder to returnFolders.
                    if any(ClassFileChanged)
                        returnFolders{retIx} = folders{ix,1};
                        returnFolders{retIx, 2} = ClassFileChanged;
                        retIx = retIx+1;
                    end
                end

                % Grow returnFolders array if necessary.
                if retIx > length(returnFolders)
                    returnFolders = [returnFolders ; cell(100,2)]; %#ok<AGROW>
                end

            end
        end

      % -- Get all folder names; remove folders starting with '.'.
        list        = dir(folders{ix,1});
        prevFolder  = folders{ix,1};
        folderIdx   = find(cellfun(@(x) x(1),{list.name}) ~= '.' & [list.isdir] == true);

        % Remove folders that cannot contain objects of current class
        pBool = false(length(folderIdx),1);
        for i = 1: length(folderIdx)
            pClassName  = regexp(list(folderIdx(i)).name,'[A-Za-z0-9]+','match','once');
            if any(strcmp(pClassName, posParents));
                pBool(i) = true;
            end
        end
        
        % Grow folders array if necessary
        if (ix2 + length(folderIdx)) > length(folders)
            folders = [folders ; cell(100,2)]; %#ok<AGROW>
        end
        
        % Add subfolders to folder array or increase ix.
        if ~isempty(folderIdx)
            folders{ix,1} = fullfile(folders{ix,1} , list(folderIdx(1)).name );
            folders{ix,2} = pBool(1);
            
            for i = 2: length(folderIdx)
                folders{ix2,1} = fullfile(prevFolder , list(folderIdx(i)).name);
                folders{ix2,2} = pBool(i);
                ix2 = ix2 + 1;
            end
        else
            ix = ix+1;
        end
        
        % Break if the active index == the first available index.
        if ix == ix2; break; end;
    end
    
    % Return the returnFolders array.
    folders = returnFolders(~cellfun(@isempty, returnFolders(:,1)),:);
    
    warning('on', 'MATLAB:class:mustReturnObject');

end
