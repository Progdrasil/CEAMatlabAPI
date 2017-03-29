classdef input_class < handle
    
    properties (SetAccess = private)
        parent; % Link to the CEA class
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