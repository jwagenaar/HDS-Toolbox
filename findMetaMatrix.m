function out = findMetaMatrix(query, pathMatrix)
    
    global HDSsearchTables

    names = {query.pattern.field};
    queryIds = zeros(length(names),1);    
    for ii = 1:length(names)
        queryIds(ii) = find(strcmp(names{ii}, HDSsearchTables.mDataTree.metaNames),1);
    end
    [queryIds loc] = sort(queryIds);
    qPatterns = {query.pattern.arg};
    qPatterns = qPatterns(loc);
    
    %Get unique query fields --> faster than matlabs version
    idxs = zeros(length(queryIds),1);
    allidx = queryIds;
    i = 1;
    while 1
        aux = allidx(find(allidx,1));
        if isempty(aux)
            break
        end
        idxs(i) = aux;
        allidx(allidx == idxs(i))=0;
        i = i +1;
    end
    UniqueQueryIds = sort(idxs(1:i-1));
    
    
    %Get unique ClassID fields --> faster than matlabs version
    idxs = zeros(length(HDSsearchTables.mDataTree.classes),1);
    allidx = pathMatrix;
    i = 1;
    while 1
        aux = allidx(find(allidx,1));
        if isempty(aux)
            break
        end
        idxs(i) = aux;
        allidx(allidx == idxs(i))=0;
        i = i +1;
    end
    classIds = sort(idxs(1:i-1));
    
    metaTable = zeros(length(classIds),length(UniqueQueryIds));
    for ii = 1: length(classIds)
        metaTable(ii,:) = ismembc(UniqueQueryIds, sort(HDSsearchTables.mDataTree.classes(classIds(ii)).metaIds));
        
    end
    
    % metaIds, classIds and queryPatterns are always sorted from low to high for use with ISMEMBC.
    
    out = struct('classIds',classIds,'metaIds',UniqueQueryIds,'table',metaTable,'qIds',queryIds,'qPatterns',[]);
    out.qPatterns = qPatterns;
   
end