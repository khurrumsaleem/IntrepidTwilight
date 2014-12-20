classdef RectangularPrism < IntrepidTwilight.AdamantWave.Shape
    
    properties
        
        % 2D principal measurements
        Length
        Width
        
        % Reference location
        Datum
        
        % 3D depth
        Depth
        
        % 2D orientation
        InclinationAngle
        
    end
    
    methods
        function RP = RectangularPrism(n)
            if (nargin >= 1) && not(isempty(n))
                RP(n) = RectangularPrism();
            end
        end
        


        function [] = SetVertices(RP,V1,V2,V3,V4)
            
            SetVertices@IntrepidTwilight.AdamantWave.Shape(RP,V1,V2,V3,V4);
            
            % Need to check for rectangular-ness
            
        end
        
        function [] = SetPrincipals(RP,Datum,Length,Width,Depth,Theta)
            
            
            RP.Length = Length;
            RP.With   = Width;
            RP.Volume = Length*Width*Depth;
            RP.Datum  = Datum;
            
        end

    end
    
end