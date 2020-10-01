classdef field_screen < handle
    %FIELD SCREEN is a class for projecting complex
    % wave fields. It can be used to represent 
    % complex field sensors such as single element photodiodes
    % for interferometry.
    % A. Schultze 01/10/2020 (GaussCAD toolbox)
    properties
        r;
        n;
        dim;
        pix;
        beams=[];
        E=[];     %complex field array
        mask=[];
        rotax = [0 0 1];
        rotang = 0;
    end
    
    methods
        function obj = field_screen(in_r,in_n,in_dim, in_pix)
            %in_r,  -[3x1] position vector
            %in_n   -[3x1] unit normal vector
            %in_dim -[2x1] dimension
            %in_pix -[2x1] pixels
            
            obj.r=in_r;
            obj.n=in_n;
            obj.dim=in_dim;
            obj.pix=in_pix;
            obj.mask = ones(in_pix);
        end
        
        function outputArg = add_beam(self,in_beam)
            % Adds a beam to list of beams
            self.beams=   [self.beams in_beam];
            outputArg = self;
        end
        
        function retval = set_mask_round(self)
            y = linspace( -self.dim(1)/2, self.dim(1)/2, self.pix(1) );
            z = linspace( -self.dim(2)/2, self.dim(2)/2, self.pix(2) );
            [ yy,zz] = meshgrid( y, z );
            ra=self.dim(1)/2;
            self.mask =((yy.^2 + zz.^2)<= ra.^2);
            retval = self.mask;
        end
        
        %sets a + shaped gap as typical for QPD diodes
        function retval = set_mask_gap(self, a_gap)
            y = linspace( -self.dim(1)/2, self.dim(1)/2, self.pix(1) );
            z = linspace( -self.dim(2)/2, self.dim(2)/2, self.pix(2) );
            [ yy,zz] = meshgrid( y, z );
            self.mask =(abs(yy)>a_gap/2 & abs(zz)>a_gap/2) & self.mask;
            retval = self.mask;
        end
        
        %gets mask which covers one quadrant
        function retval = mask_quadrant(self,a_q)
            y = linspace( -self.dim(1)/2, self.dim(1)/2, self.pix(1) );
            z = linspace( -self.dim(2)/2, self.dim(2)/2, self.pix(2) );
            [ yy,zz] = meshgrid( y, z );
            
            %      4 | 1      ^(y)
            %     -----       | 
            %      3 | 2      --->(z)
            switch (a_q)
                case 1
                    qq = (yy>0) & (zz > 0); %+y +z
                case 2
                    qq = (yy<0) & (zz > 0); %-y +z
                case 3
                    qq = (yy>0) & (zz < 0); %+y -z
                case 4
                    qq = (yy<0) & (zz < 0); %-y -z
                otherwise
                    warning('No valid quadrant specified');
            end 
            retval = qq;
        end
        
        
        function [int,ph,e] = calc_int_phase(self)
            int=abs(sum(self.E.*self.mask,'all')); %total intensity
            ph =angle(sum(self.E.*self.mask,'all')); %total phase angle
            e  = sum(self.E.*self.mask,'all'); %complex field response
        end
        
        function [int,ph,e] = calc_int_phase_quadrants(self)
            int=zeros(1,4);
            ph=zeros(1,4);
            e =zeros(1,4);
            for i=1:4
                maskq = self.mask_quadrant(i);
                int(i)=abs(sum(self.E.*(maskq & self.mask),'all')); %total intensity
                ph(i) =angle(sum(self.E.*(maskq & self.mask),'all')); %total phase angle
                e(i)  = sum(self.E.*(maskq & self.mask),'all');
            end
        end
        
        function [dws_y,dws_z] = calc_dws(self)
            %TODO: Check definition of y/z
            [~,ph,e] = calc_int_phase_quadrants(self);
            dws_y =  (ph(1)+ph(4))-(ph(2)+ph(3));
            dws_z =  (ph(1)+ph(2))-(ph(3)+ph(4));
            %dws_y = angle(e(1)+e(4))-(e(2)+e(3));
            %dws_z = angle(e(1)+e(2))-(e(3)+e(4));
            
        end
        
        function [ac,dc,phi,int] = calc_contrast(self)
            %figure out ac and dc of interference
            %by varying phase of one beam and find
            %min and max
            
            ni=36;
            for i=1:ni
                phi(i)=2*pi/i;
                self.shift_phase_beam1( phi(i));
                int(i)=sum(abs(sum(self.E,3).^2).*self.mask,'all');
            end   
            
            ac = max(int)-min(int);
            dc = min(int);
        end
        
        function retval = plot_mask(self)
            y = linspace( -self.dim(1)/2, self.dim(1)/2, self.pix(1) );
            z = linspace( -self.dim(2)/2, self.dim(2)/2, self.pix(2) );
            [ yy,zz] = meshgrid( y, z );
            subplot(2,2,1);
            retval=surf(yy,zz,double(self.mask),'EdgeColor','none');   
            xlabel('y');ylabel('z');
            view([0 90]);
            for i=1:3
               subplot(2,2,1+i); 
               surf(yy,zz,double(self.mask & self.mask_quadrant(i)),'EdgeColor','none');
               xlabel('y');ylabel('z');
               legend(['Q' num2str(i)],'Location','southeast');
               view([0 90]);
            end      
        end
        
        function retval = render(self)
            % Pixels and their position in space
            % Without rotation Reference Y = Screen Z ("up");

            y = linspace( -self.dim(1)/2, self.dim(1)/2, self.pix(1) );
            z = linspace( -self.dim(2)/2, self.dim(2)/2, self.pix(2) );
            x = 0;
            [ xx,yy, zz ] = meshgrid( x, y, z );
            %points = [xx(:), yy(:) zz(:)];
            
            S = [ xx(:) yy(:) zz(:) ];
            % rotate and shift
            if self.rotang ~= 0
                S = rodrigues_rot( S, self.rotax, self.rotang );
            end     
            %calculate the E field for each position for each beam           
            for i=1:length(self.beams)
                this_beam=self.beams(i);
                this_E=this_beam.calc_field_xyz( S + self.r);
                self.E(:,:,i)= reshape(this_E,size(xx));
            end
            retval= self.E;                      
        end % function retval = render(self)
        
        %This function will simplify and speed up
        %shifting phase of one beams rendered result
        function ret_handle = shift_phase_beam1(self,a_phi)
            self.E(:,:,1)=exp(-1i*a_phi).*self.E(:,:,1);    
        end
        
        function ret_handle = plot_intensity(self)
            y = linspace( -self.dim(1)/2, self.dim(1)/2, self.pix(1) );
            z = linspace( -self.dim(2)/2, self.dim(2)/2, self.pix(2) );
            [ xx,yy ] = meshgrid( y,z );
            subplot(2,2,1);
            ret_handle=surf(xx,yy,sum(abs(self.E),3),'EdgeColor','none');
            title('Sum of Intensity');
            subplot(2,2,3);
            
            for i=1:size(self.E,3)
                ret_handle=surf(xx,yy,abs(self.E(:,:,i)),'EdgeColor','none');hold on;
            end
            title('Each Intensity');
            legend();
            
            figure();
            if size(self.E,3) == 2
                ac = min(abs(self.E),[],3);
                dc = sum(abs(self.E),3)-ac;
                subplot(2,1,1);
                surf(xx,yy,ac,'EdgeColor','none');hold on;
                legend('interfering');
                subplot(2,1,2);
                surf(xx,yy,dc,'EdgeColor','none');hold on;
                legend('non-interfering');
            end
        end %plot intensity
        
        
        function ret_handle = plot_interference(self)
            y = linspace( -self.dim(1)/2, self.dim(1)/2, self.pix(1) );
            z = linspace( -self.dim(2)/2, self.dim(2)/2, self.pix(2) );
            [ xx,yy ] = meshgrid( y,z );
            subplot(2,2,1);
            %int=abs(sum(self.E,'all')); %total intensity
            %ph =angle(sum(self.E,'all')); %total phase angle
            
            %max possible intensity for scaling
            max_int = max(sum(abs(self.E),3).^2,[],'all'); 
            
            subplot(2,1,1)
            ret_handle=surf(xx,yy,abs(sum(self.E,3)).^2,'EdgeColor','none');
            zlabel('Intensity');
            colorbar; view([0,90]); caxis([0 max_int]);
            xlabel('Sensor Y (mm)');ylabel('Sensor Z (mm)');
            title('Interference (Intensity)');
            subplot(2,1,2)
            ret_handle=surf(xx,yy,abs(angle(sum(self.E,3))),'EdgeColor','none');
            zlabel('Phase (rad)');
            title('Interference (Phase)');
            colorbar;caxis([0 pi]);view([0,90]);
            xlabel('Sensor Y (mm)');ylabel('Sensor Z (mm)');
        end % plot interference
        
    end
end

