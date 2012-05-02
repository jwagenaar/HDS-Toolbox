function out = formQueryStruct(returnclass, treeIx, queryCell)
    %FORMQUERYSTRUCT  Formats the query in a standard form the search
    %funtion accepts.
    %
    
    % Comparison values are defined can be:
    % '==' '~=' '<=' '<' '>' '>=' '<>' and'[]' where '<>' means 'to include any'
    % and '[]' means to include all. This applies to numbers or strings.

    global HDSsearchTables
    global HDSManagedData


    % Iterate over queryCell and extract query structure
    ix = 1; %index in the queryCell
    jx = 0; % index in queryStruct
    expr = '(?<field>[A-Za-z]+)\s+(?<comp>[<>[]=~]+)\s+(?<arg>[A-Za-z .0-9]*)';
    qStruct = struct('field',[],'comp',[],'arg',cell(10,1));
    while ix <= length(queryCell)
        
        jx = jx+1;

        % First index should contain a property and a comparator.
        try
            qStruct(jx) = regexp(queryCell{ix}, expr, 'names');
        catch ME
            throwAsCaller(MException('HDS:HDSSFIND','HDSFIND: Unable to parse query syntax.'));
        end
            
            
        % If ix is not last index of string, find argument in sam string
        if isempty(qStruct(jx).arg)
            if length(queryCell) >= ix+1;
                qStruct(jx).arg = queryCell{ix+1};
                ix = ix+2;
            else
                throwAsCaller(MException('HDS:HDSFIND','HDSFIND: Incorrect query syntax.'));
            end
        else
            if all(double(qStruct(jx).arg) >= 46 & double(qStruct(jx).arg) <= 57)
                qStruct(jx).arg = str2double(qStruct(jx).arg);
            elseif strcmp(strtrim(qStruct(jx).arg),'true')
                qStruct(jx).arg = true;
            elseif strcmp(strtrim(qStruct(jx).arg),'false')
                qStruct(jx).arg = false;
            end
            ix = ix+1;
        end
        
    end
    
    
    % Load the HDSsearchTables if not previously loaded.
    if isempty(HDSsearchTables)
        loadsearchtables(HDSManagedData(treeIx).basePath, HDSManagedData(treeIx).treeConst(1));
    elseif HDSManagedData(treeIx).treeConst(1) ~= HDSsearchTables.hostId
        loadsearchtables(HDSManagedData(treeIx).basePath, HDSManagedData(treeIx).treeConst(1));
    end
  
    out = struct('return',returnclass, 'pattern', qStruct(1:jx), 'metaClasses',[], 'metaClassIdx',[]);

    % AllClasses is boolean array that identified all classes that contain
    % one or more propertynames specified in the query
    allClasses = false(length(HDSsearchTables.mDataTree.classes),1);
    
    
    % FieldIds are the indeces of the metanames in the
    % HDSsearchTables.mDataTree.metaNames
    fieldIds = zeros(length(out.pattern),1);
    
    for ii = 1:length(fieldIds)
        if strcmp('link',out.pattern(ii).field)
            % Trying to match to linked objects. 
            out.pattern(ii).field = out.pattern(ii).field;
            
            aux = [out.qStruct(ii).arg.objIds];
            aux = [aux{:}];
            out.pattern(ii).value = aux(1,:);
            out.pattern(ii).comp = '<>';
            fieldIdx = 0;
        else
            
            fieldIdx = find( strcmp( out.pattern(ii).field, HDSsearchTables.mDataTree.metaNames),1);

            if isempty(fieldIdx)
                throwAsCaller(MException('formQueryStruct:UnkownField',sprintf('Unknown metafield: %s.',out.pattern(ii).field)));
            end
            
        end
        
        
        
        fieldIds(ii) = fieldIdx;
    end

    for iClass = 1: length(allClasses)
        allClasses(iClass) = any(ismembc(fieldIds, HDSsearchTables.mDataTree.classes(iClass).metaIds));
    end
    
    % Out.metaClasses is a list of all classes that have at least one of
    % the meta search terms are located in the specified class. This is
    % used to get all possible search paths. 
    out.metaClasses = {HDSsearchTables.mDataTree.classes(allClasses).name};
    out.metaClassIdx = find(allClasses);
    
end