classdef input_class < handle
    
    properties (SetAccess = private)
        parent;
    end
    
    methods (Access = ?CEA)
        % Constructeur
        function this = input_class(parent)
            this.parent = parent;
        end
    end
    methods
        % Autres fonctions
        inp = rocket(obj)
    end
end