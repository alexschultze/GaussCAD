classdef bench_abcd < handle
    %This class groups element of class "element" which describe
    %abcd matrices.
    %By using the absolute distance a complete system is built
    %by adding free space distances
    
    properties
        elements=[];
        elements_pos=[];
        lambda =532e-9; %(m) wavelength
    end
    
    methods (Static)
        function ret_q=propagate_gauss(abcd,q)
            % 1,1 - A
            % 1,2 - B
            ret_q=(abcd(1,1).*q+abcd(1,2))./ (abcd(2,1).*q+abcd(2,2));
        end

        function [R]=gauss_curvature(z,zr)
            R=z.*(1+(zr./z).^2);  
        end
        
        function [R]=gauss_curvature_q(q)
            R=real(q).*(1+(imag(q)./real(q)).^2);  
        end

        
        
    end
    
    methods

        function obj = bench_abcd(par_lambda)
            % Construct an instance of this class
            %   Detailed explanation goes here
             if exist('par_lambda','var')
                obj.lambda = par_lambda;
             end
        end

        function [w0]=calc_w0(obj,q)
            w0=sqrt(imag(q)*obj.lambda/pi);
        end       

        function [w]=calc_w(obj,q)
            w = sqrt(imag(q)*obj.lambda/pi).*sqrt(1+(real(q)./imag(q)).^2);
        end


        function add(obj, in_pos,in_element)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.elements=[obj.elements in_element];
            obj.elements_pos    =[obj.elements_pos in_pos];
        end
        
        function plot(obj, par_rl)
            
             if(size(par_rl)~=[2 1])
                 error('Wrong parameter dimension [r alpha]');
             end
             %figure();
             subplot(2,1,1);
             hold on;
             %draw elements
             
             
            for i=1:size(obj.elements,2)
                xline(obj.elements_pos(i),'--r');
                
                if(i>1)
                    last_z=obj.elements_pos(i-1);
                    last_rl=this_rl(:,i-1);
                else
                    last_z=0;
                    last_rl=par_rl;
                end

                z= obj.elements_pos(i)-last_z;
                abcd = [1,z;0,1];
                tmp_rl=abcd*last_rl; %free space
                this_rl(:,i)=obj.elements(i).abcd*tmp_rl; 
                

                
            end
            
            % Build complete with initial
            rl= [par_rl, this_rl];
            pos=[0,obj.elements_pos];
            plot(pos, rl(1,:));
            xlabel('z (m)');
            ylabel('r(z) (m)');
            subplot(2,1,2);
            stairs(pos, rl(2,:));
            xlabel('z (m)');
            ylabel('alpha(z) (rad)');
            sgtitle('Matrix Optics: Paraxial Beam Propagation');
        end
 

        
        function [q,pos, R,w,w0]=plot_gauss(obj, in_q)

             %draw elements
             
             
            for i=1:size(obj.elements,2)              
                if(i>1)
                    last_z=obj.elements_pos(i-1);
                    last_q=this_q(:,i-1);
                else
                    last_z=0;
                    last_q=in_q;
                end

                z= obj.elements_pos(i)-last_z;
                abcd = [1,z;0,1];
                tmp_q=obj.propagate_gauss(abcd,last_q);
                this_q(i)=obj.propagate_gauss(obj.elements(i).abcd,tmp_q);
                
                 %add two elements for each calculation, before
                %and after element
                pos(2*i-1)=obj.elements_pos(i)-eps;
                q_all(2*i-1)=tmp_q;
                pos(2*i)=obj.elements_pos(i);
                q_all(2*i)=this_q(i);       
                
                %Radius of curvature R=z(1+(zr/z).^2)
                R(2*i-1) =obj.gauss_curvature_q(tmp_q);
                R(2*i)   =obj.gauss_curvature_q(this_q(i));
                
                
            end
            
            
            % Build complete with initial
                q= [in_q, q_all];
                pos=[0,pos];
                R =[obj.gauss_curvature_q(in_q), R];
            % Calculate corresponding beam widths
                w = obj.calc_w(q);
                w0 = obj.calc_w0(q);
            if(nargout==0)
                subplot(4,1,1);
                hold on;

                %draw each element
                for i=1:length(obj.elements_pos)
                    xline(obj.elements_pos(i),'--r');
                    if ~strcmp(obj.elements(i).type,'translation')
                        text(obj.elements_pos(i),0,obj.elements(i).type,'Rotation',90,'Color','r','Interpreter','none');
                    end
                
                end
                xlim([0,obj.elements_pos(end)]);

                xlabel('Setup (m)');

                subplot(4,1,4);
                stairs(pos, imag(q));
                xlabel('z (m)');
                ylabel('z0(z) (rad)');
                
                subplot(4,1,3);
                plot(pos, R,'*--');
                xlabel('z (m)');
                ylabel('R(z) (m)');

                subplot(4,1,2);
                plot(pos, w,'--');
                xlabel('z (m)');
                ylabel('w(z) (m)');

                sgtitle('Matrix Optics: Gauss Beam Propagation');
                %replace suptitle in >=r2018b
            end
        end 
        
    end
end

