classdef HDSFileRepos < handle

  properties (SetAccess = private)
    typeID   = '' % Type of the repository, restricted options
    typeAttr = {} % Attributes for type, depending on type definition.
    root = ''     % Root of repository (usually folder of the files)
    chNames = {}  % Channel Names.
    sf = 0        % Sampling frequency.
    files = {}    % FileNames Rows are channels, columns are blocks
  end

  properties (Transient)
    userData    = {}  % Can be used by getMethod to store stuff in object
    fetchCache  = []  % Holds data if necessary.
  end
  
  methods
    function obj = HDSFileRepos(root, chNames, sf, files, type, varargin)
      assert(any(strcmp(type,{'BinByChannel' 'UINT32BINBYBLOCK'})), ...
        'Uncompatible Type');
      obj.typeID = type;
      obj.root = root;
      obj.chNames = chNames;
      obj.sf = sf;
      obj.files = files;
      
      if nargin > 5
        names = varargin(1:2:(end-1));
        values = varargin(2:2:end);
        for i = 1: length(names)
          obj.typeAttr.(names{i}) = values{i};
        end
      end
      
        
     
    end

    function data = getData(obj, channels, indeces)
      switch obj.typeID
        case 'BinByChannel'
          data = getBinByChannel(obj, channels, indeces);
          
      end
    end
    
    function attr = getAttr(obj)
      
      attr = struct('chNames',[], 'sf' , obj.sf);
      attr.chNames = obj.chNames;
      switch obj.typeID
        case 'BinByChannel'
          attr = attrBinByChannel(obj, attr);
          
      end
    end
    
    data = getuintbinbychannel(obj, channels, indeces);
    
  end
    
  
end