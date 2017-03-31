classdef input_class < handle
    % input_class Creates ioinp with CEA parameters
    %   Gets parent CEA object and its parameters and creates the cell array
    %   of strings which will be the input for the mex cea function. it will
    %   be created for which ever method is called (application).
    %   Currently the only applications supported in the mex file is rockets.
    % 
    % input_class properties:
    %   parent - The parent CEA class for this input
    %
    % input_class methods:
    %   input_class - constructor that only the CEA class can access
    %   rocket - This method creates ioinp for the Rocket application
    %
    % input_class methods not yet supported:
    %   hp
    %   tp
    %   det
    %   sp
    %   tv
    %   uv
    %   sv
    %   shock
    %
    % See Also:
    %   CEA
    properties (SetAccess = private)
        % Link to the parent CEA class
        % 
        % See Also:
        % CEA
        parent; % Link to the parent CEA class
    end
    
    methods (Access = ?CEA)
        % Constructeur
        function this = input_class(parent)
            % Creates link to parent CEA class object
            this.parent = parent; 
        end
    end
    methods
        % Autres fonctions
        inp = rocket(obj)
    end
end