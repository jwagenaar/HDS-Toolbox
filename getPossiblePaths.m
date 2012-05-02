function out = getPossiblePaths(metaMatrix, pathMatrix)
    %GETPOSSIBLEPATHS  Checks pathmatrix for possible paths
    %
    %   GETPOSSIBLEPATHS

    metaMatrixLength =  length(metaMatrix.metaIds);
    pathMetaTable = zeros(size(pathMatrix,1), metaMatrixLength);
    pathMatrixSize = size(pathMatrix);
        
    for ii = 1:size(pathMatrix,1)
        idx = 1;
        while idx <= pathMatrixSize(2)
            checkClassIx = pathMatrix(ii,idx);
            if ~checkClassIx
                break
            end
            pathMetaTable(ii,:) = pathMetaTable(ii,:) | metaMatrix.table(metaMatrix.classIds == checkClassIx ,:);
            idx = idx + 1;
        end
    end
    
    out = find(all(pathMetaTable,2));
end