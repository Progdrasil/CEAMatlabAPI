classdef CEA < handle

    properties (SetAccess = public)
        data; % Struct output of the mex function
        ioinp = 'No inputs created'; % Input cell array of strings
        % Input_Class object to generate ioinp
        % See Also input_class.m
        input; % Input_Class object to generate ioinp

        %Parameters for creation of input
        
        OF; % Oxidizer to fuel ratio
        pressure; % Input pressure
        presUnit = 'psia'; % Pressure unit
        supar; % Super sonic area ratio
        PcPe; % Pc/Pe
        Debug = false; % to have class debug outputs

    end
    
    properties (Access = ?input_class)
        fuel; % Fuel names cell vector of strings
        fuelWt; % Fuel weight ratio vector [%%]
        fuelTemp; % Fuel tempuratures vector

        oxid; % Oxidizer names cell vector of strings
        oxidWt; % Oxidizer weight ratio vector [%%]
        oxidTemp; % Oxidizer tempuratures vector
    end
    methods
        % Constructeur
        function this = CEA()
            % calls input_class constructor for the input property
            this.input = input_class(this);
        end
        
        data = run(obj) % Runs the compiled mex function
        setFuel(obj, fuels, fuelWeights, fuelTemps) % Function to set fuel property's
        setOxid(obj, oxids, oxidWeights, oxidTemps) % Function to set oxidizer property's
    end
end