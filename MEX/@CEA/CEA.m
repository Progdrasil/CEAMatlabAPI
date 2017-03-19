classdef CEA < handle

    properties (SetAccess = public)
        data;
        ioinp = 'No inputs created';
        input;
        % variables for Input creation
        OF;
        pressure;
        presUnit = 'psia';
        supar;
        PcPe;
        Debug = false;

    end
    
    properties (Access = ?input_class)
        fuel;
        fuelWt;
        fuelTemp;

        oxid;
        oxidWt;
        oxidTemp;
    end
    methods
        % Constructeur
        function this = CEA()
            this.input = input_class(this);
        end
        
        % Autres fonctions
        data = run(obj)
        setFuel(obj, fuels, fuelWeights, fuelTemps)
        setOxid(obj, oxids, oxidWeights, oxidTemps)
    end
    
end