function loadsearchtables(path, hostId)
    global HDSsearchTables

    HDSsearchTables = struct('mDataTree','','hostId',hostId);
    
    HDSsearchTables.mDataTree = load(fullfile(path,'HDSsearch','HDS_Treestruct.mat'));
    for iClass = 1 : length(HDSsearchTables.mDataTree.classes)
        HDSsearchTables.(HDSsearchTables.mDataTree.classes(iClass).name) = load(fullfile(path,'HDSsearch',sprintf('S%s.mat',HDSsearchTables.mDataTree.classes(iClass).name)));
    end
    
    
end
