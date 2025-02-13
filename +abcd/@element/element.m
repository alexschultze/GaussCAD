classdef element < handle
    %ELEMENT is a Paraxial element
    % described by ABCD matrix
    % A. Schultze 01/10/2020 (GaussCAD toolbox)
    properties
        abcd;
        type;
    end
    
    methods
        
        function obj = disp(obj)
            disp(obj.type);
            disp('[ABCD]: ');
            disp(obj.abcd);
        end
        
        function obj = element(in_type,varargin)
            %Construct an instance of this class
            obj.type=in_type;
            switch lower(char(in_type))
                case 'translation'
                    d=varargin{1};
                    obj.abcd = [1,d;0,1];
                case 'lens'
                    f=varargin{1};
                    obj.abcd = [1,0;1/-f,1]; 
                case 'lens_thick'
                    n1=varargin{1};
                    n2=varargin{2};
                    R1=varargin{3};
                    R2=varargin{4};
                    
                    %obj.abcd = [1,0;1./-f,1]; 
                case 'refraction'
                    n1=varargin{1};
                    n2=varargin{2};
                    obj.abcd = [1,0;0,n1./n2]; 
                case 'mirror'
                    obj.abcd = [1,0;0,1]; 
                case 'mirror_curved'
                    p=1/varargin{1};
                    obj.abcd = [1,0;-2*p,1];       
                case 'spherical_refraction'
                    n1=varargin{1};
                    n2=varargin{2};
                    p=varargin{3};
                    obj.abcd = [1,0;(n1./n2-1).*p,n1./n2]; 
                case 'screen'
                    obj.abcd = [1,0;0,1]; 
                otherwise
                   obj.abcd=varargin{1}; 
            end %switch
            
        end %fcn
        
        %Propagate the ray vector (r,phi)
        function ret_r = propagate(self,in_r)
            %Calculate new r matrix
            ret_r = self.abcd*in_r;
        end
        
        %Propagate the complex beam parameter (Gauss)
        function ret_q = propagate_gauss(self,in_q)
        ret_q=(self.abcd(1,1).*in_q+self.abcd(1,2))./ (self.abcd(2,1).*in_q+self.abcd(2,2));
        end
        
    end
end

