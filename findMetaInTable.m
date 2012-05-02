function resultMatrixOut = findMetaInTable(mHDS, classId, resultMatrixIn, metaMatrix, path, elim)
    %FINDMETAINTABLE returns indeces in class that satisfy query.
    %
    %   RESULTOUT = FINDMETAINTABLE(mHDS, CLASSID, RESULTIN, METAMATRIX,
    %   PATH, ELIM) 
    %
    %   If ELIM is true, the metaMatrix is checked whether the resulting
    %   indeces can still result in a match after each iteration. This
    %   speeds up large searches when 
    %
    %   Always search down the tree to eliminate as many objects as
    %   possible early in the iterations.

    % Get info on CLASS that is being searched
    qTableRow       = ismembc2(metaMatrix.qIds, metaMatrix.metaIds);        % Indicates which rows in table correspond with which metaId
    qFieldNames     = mHDS.mDataTree.metaNames(metaMatrix.qIds);            % Indicated the names of the fields of query.
        
    mFieldRow       = find(metaMatrix.table(metaMatrix.classIds == classId,:)); % Which metaIds are in the current class? 
    fieldInClass    = find(ismembc(qTableRow, mFieldRow));                      % Which colums in table are in the current class?
    
    %Get CountRow
    countRows = zeros(length(metaMatrix.metaIds),1);
    for ii = 1:length(fieldInClass)
        countRows(qTableRow(fieldInClass(ii))) = countRows(qTableRow(fieldInClass(ii))) + 1;
    end
    
    % Get correct Table
    table = mHDS.(mHDS.mDataTree.classes(classId).name);
    
    % Find index in path
    pathIdx = find(classId==path,1);
    
    if isempty(resultMatrixIn)
       %resultMatrixIn = struct('classIdx',1:length(table),'matrix',false(length(table), length(metaMatrix.qIds)));
       resultMatrixIn = struct('classIdx',1:length(table.parentId),'matrix',false(length(table.parentId), length(metaMatrix.qIds)));
    end
    
    
    % -- Prepare search by eliminating all objects with incorrect parent--

    % Check resultMatrixIn
    if size(resultMatrixIn.matrix, 2) ~= length(metaMatrix.qIds)
        throw(MException('FindMetaInTable:inCorrect_input','The RESULTMATRIXIN has incorrect dimensionss.'));
    end

    % Find order in which to find meta
    includeParentIdxI = true(length(resultMatrixIn.classIdx),1);
    if pathIdx ~= length(path)
        tempTable = metaMatrix.table(ismembc(metaMatrix.classIds, sort(path(pathIdx +1:end))), :);    

        if isempty(tempTable)
            metaCount = zeros(1,size(tempTable,2));
        else
            metaCount = sum(tempTable,1);
        end

        for ii = 1:length(metaCount)
            metaCount(ii) = metaCount(ii) + countRows(ii);
            if metaCount(ii) == 0
                %No possibility to meet requirement for this metatag in
                %query.
                includeParentIdxI = includeParentIdxI & resultMatrixIn.matrix(:,ii);
            end
            resultMatrixIn.matrix = resultMatrixIn.matrix(includeParentIdxI,:);
            resultMatrixIn.classIdx = resultMatrixIn.classIdx(includeParentIdxI);
        end
    else
        metaCount = zeros(length(qTableRow));
        for ii = 1:length(metaCount)
            metaCount(ii) = metaCount(ii) + countRows(ii);
        end

    end

    % Elim by classId & ClassIdx
    if elim
        comboIds = table.comboId;
        parentIds = single(resultMatrixIn.classId *1e-4 + resultMatrixIn.classIdx);

        includeParentIdx    = ismembc2(comboIds, sort(parentIds) );  % Find ids of parents and logical of combo
        includeIdx          = logical(includeParentIdx);

        allidx = includeParentIdx;
        includeParentIdx = unique(allidx);

        includeParentIdx(includeParentIdx==0)=[];


        % Build result Matrix
        result = false(sum(includeIdx), length(metaMatrix.qIds));

        aux = comboIds(includeIdx);
        aux2 = resultMatrixIn.matrix( includeParentIdx,:);
        for ii = 1: length(includeParentIdx)
            ix = aux == parentIds(includeParentIdx(ii));            
            result(ix,:) = aux2(ii*ones(sum(ix),1), :);
        end

    else
        %includeIdx = true(1,length(table));
        includeIdx = true(1,length(table.parentId));
        result = false(sum(includeIdx), length(metaMatrix.qIds));
    end


    

    %  --- Compare fields in class to the query input ---
    activeResultRows = true(size(result,1),1);
    for ii = 1: length(fieldInClass)
        
        if elim
            includeIdxII  = includeIdx;
            includeIdxIII = activeResultRows & ~result(:, fieldInClass(ii));
            includeIdxII(includeIdxII) = includeIdxIII;
        else
            includeIdxII  = includeIdx;
            includeIdxIII = 1:size(result,1);
        end
        
        if any(includeIdxIII)
        
            metaFields = table.(qFieldNames{fieldInClass(ii)});
            metaFields = metaFields(:,includeIdxII);

            pattern = metaMatrix.qPatterns{fieldInClass(ii)};

            if ischar(metaFields)
                aux = result(includeIdxIII, fieldInClass(ii));
                if ischar(pattern)
                    sizeDiff = size(metaFields,1) - length(pattern);
                    if sizeDiff >=0
                       sPat = pattern(ones(size(metaFields,2),1),:)';
                       a = ' ';
                       sPat = [sPat ; a(ones(sizeDiff, size(metaFields,2)))]; %#ok<AGROW>
                       tRes = ~any(sPat - metaFields,1);
                       result(includeIdxIII, fieldInClass(ii)) = aux | tRes';
                        
                    end
                else
                    throwAsCaller(MEXception('HDS:hdsfind',sprintf('HDSFIND: Incorrect class search pattern; pattern is of class: %s, values are of class: CHAR.',upper(class(pattern)))));
                end
                
            elseif isnumeric(pattern) || islogical(pattern)  
                aux = result(includeIdxIII, fieldInClass(ii));
                result(includeIdxIII, fieldInClass(ii)) = aux | pattern == metaFields' ;  
                
            elseif iscell(metaFields)
                res = false(length(metaFields),1);
                for jj=1: length(metaFields)
                    
                    if iscellstr(metaFields{jj})
                        res(jj) = any(strcmp(pattern, metaFields{jj}));
                    end
                end
                aux = result(includeIdxIII, fieldInClass(ii));

                result(includeIdxIII, fieldInClass(ii)) = aux | res ;

            end

            metaCount(qTableRow(fieldInClass(ii))) = metaCount(qTableRow(fieldInClass(ii)))-1;
            if metaCount(qTableRow(fieldInClass(ii))) == 0 
                activeResultRows = activeResultRows & result(:, fieldInClass(ii));
            end
            
        end
        
    end
    
    includeIdx(includeIdx) = activeResultRows;
    pId = table.comboId(includeIdx);
    classIdx = find(includeIdx);
    cId = single(classId)*1e-4 + single(classIdx);

    resultMatrixOut = struct('classId',classId,'classIdx',classIdx,'cId',cId,'pId',pId,'matrix',result(activeResultRows,:));
        

end