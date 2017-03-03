classdef CEA

    properties (SetAccess = private)
        data
        input

    end
    methods (MethodAttributes)
        function data = run(obj)
            [pathstr,~,~] = fileparts(currentpath);
            obj.data = dos(strcat(pathstr,'/PCEA2.out'));
            
        end
        function data = get.data(obj);
            
        end
        

    end
    
end