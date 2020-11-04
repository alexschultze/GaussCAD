classdef gauss_beam < handle & matlab.mixin.Copyable
    %This class defines a gauss beam
    % A. Schultze 01/10/2020 (GaussCAD toolbox)
    properties
        
        % A gauss beam can be fully described by 3 independent parameters.
        % zr = nrefr*pi*wo.^2 / lambda0
        
        zr;         % Rayleigh length
        lambda;     % Wavelength
        w0;         % Beam Waist
        nrefr;      % Refractive Index
        
        p=[];       % Position Vector
        n=[];       % Direction Vector
        
        E0=1;       %([])  initial intensity
        phase0 = 0; %(rad) initial phase
    end
    
    methods (Static)
        function o_zr = calc_zr(a_wo, a_n,a_lambda)
            % {\displaystyle z_{\mathrm {R} }={\frac {\pi w_{0}^{2}n}{\lambda }}.}
            o_zr = pi*a_wo.^2*a_n ./ a_lambda;
        end
        
         function o_w0 = calc_w0(a_zr, a_n,a_lambda)
            o_w0 = sqrt(a_zr .* a_lambda ./( a_n .* pi));
         end
    end
    
    
    methods
            % This is the default constructor which takes 3 parameters
            % and calculates the forth (z_r).
            % Alteratively: If a forth parameter is added, 
            % w0 is ignored and calculated  from zr.
            % This constructor [p] defines the position of the beam waist!
        function obj = gauss_beam(a_pv, a_nv, a_w0,a_n,a_lambda,varargin)
            % 4 parameters: w0 / zr / k / lambda
            obj.lambda = a_lambda;
            obj.nrefr = a_n;
                        
            if length(varargin)==0
                obj.w0 = a_w0;
                obj.zr = gaussian.gauss_beam.calc_zr(a_w0,a_n,a_lambda); 
            else
                obj.zr=varargin{1};
                obj.w0 = gaussian.gauss_beam.calc_w0(obj.zr,a_n,a_lambda); 
            end
            
            obj.p = a_pv;
            obj.n = a_nv;
            obj.E0=1;      
        end
        
        function retval = calc_field_rz(self,z,r)  
            [~,indx]=find(z==0);
             z(indx)=eps;
             
            wz = self.w0.*sqrt(1+(z./self.zr).^2);
            Rz = z.*(1+(self.zr./z).^2 );
            gouy = atan(z./self.zr);
            k= 2*pi*self.nrefr./self.lambda;
            
            E = self.E0 * self.w0./wz .* ...
                exp(-r.^2./wz.^2).*...
                exp(-1i.*(k.*z+k.*r.^2./(2*Rz)- gouy));
            E = E.*exp(1i*self.phase0);
            retval = E;
        end

        function retval = calc_field_xyz(self,p)
            
                [ m, n ] = size( p );
                if n ~= 3
                    error( 'Input vector is not 3 dimensional!' );
                end
            
               cnt=size(p,1);
               %transform x,y,z into r,z
               % orthogonal distance z
               do=dot(repmat(self.n,cnt,1)',(p-self.p)')./norm(self.n);             
               % perpendicular distance r
               dp=vecnorm(cross(repmat(self.n,cnt,1),(p-self.p)),2,2 )./norm(self.n );
               retval=self.calc_field_rz(do,dp');        
        end
        
        %This will get the complex beam parameter q
        %Because it is based on paraxial optics
        %It needs a distance on axis from a given point
        function ret_q = calc_q(self,zd)
           %zd = (z-z0), current distance from beam waist
           % (-) before beam waist (+) after beam waist
           ret_q=(zd)+1i*self.zr;
        end

            %this transforms the beam as if the optical element is
            %at the origin [0 0 0]
        function ret_beam = transform_beam(self,a_q_old, a_q_new)
            ret_beam = copy(self); %copy
            
            %beam parameter
            ret_beam.zr = imag(a_q_new);
            ret_beam.w0 = self.calc_w0(ret_beam.zr,self.nrefr, self.lambda);
            
            zdiff = real(a_q_new-a_q_old); %position of beam waist old vs new
            ret_beam.p = ret_beam.p+self.n * zdiff;
        end
        
        
    end
end

