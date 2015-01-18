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
        
        function [] = SetPrincipals(RP,Datum,Length,Width,Depth,Theta)
            
            RP.Length = Length;
            RP.Width  = Width;
            RP.Volume = Length*Width*Depth;
            RP.Datum  = Datum;
            
            Ones = ones(size(Datum,1));
            Zeros = Ones * 0;
            
            %   Canonical vertex placement
            V1 = [Zeros       , -Width/2 * Ones];
            V2 = [Length*Ones , -Width/2 * Ones];
            V3 = [Length*Ones , +Width/2 * Ones];
            V4 = [Zeros       , +Width/2 * Ones];
            
            %   Setup rotation matrix
            CosTheta = cos(Theta);
            SinTheta = sin(Theta);
            Rotate   = @(V) ...
                [   CosTheta.*V(:,1) - SinTheta.*V(:,2),...
                    SinTheta.*V(:,1) + CosTheta.*V(:,2)];

            %   Translate and rotate
            V1 = Datum + Rotate(V1);
            V2 = Datum + Rotate(V2);
            V3 = Datum + Rotate(V3);
            V4 = Datum + Rotate(V4);
            
            % Put vertices into the (Vx,Vy) shape couple with ccw
            % orientation
            RP.Vx = [V1(:,1),V2(:,1),V3(:,1),V4(:,1)];
            RP.Vy = [V1(:,2),V2(:,2),V3(:,2),V4(:,2)];
            
        end
        
        function h = plot(RP,varargin)
            h = plot(RP.Vx([1:end,1]),RP.Vy([1:end,1]),varargin{:});
        end

    end
    
end