classdef CEA < handle
    % CEA Summary
    % Class to store inputs and outputs of NASA's CEA program and run it
    %
    % CEA Properties:
    %   OF - Oxidizer to fuel ratio
    %   pressure - Input chamber pressure
    %   presUnit - Pressure Unit psia by default
    %   supar - Supersonic area ratio
    %   PcPe - Pc/Pe
    %   input - input_class object to genereate ioinp
    %   ioinp - Input cell array of strings passed to the cea mex function 
    %   data - Struct output of the cea mex function
    %
    % CEA Private properties:
    % Set by using the setFuel and setOxid methods
    %   fuel - Fuel names stored in a cell vector of strings
    %   fuelWt - Mass weight fraction of the fuels in a vector
    %   fuelTemp - Temperature of fuels in a vector
    %   oxid - Oxidizer names stored in a cell vector of strings
    %   oxidWt - Mass weight fraction of the oxidizers in a vector
    %   oxidTemp - Temperature of oxidizers in a vector
    %
    % CEA Methods:
    %   setFuel - Set the parameters fuel, fuelWt and fuelTemp correctly
    %   setOxid - Set the parameters oxid, oxidWt and oxidTemp correctly
    %   run - Run's the cea mex function with the current ioinp parameter

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
        supar; % Supersonic area ratio
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