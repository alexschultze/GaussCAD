classdef bench_abcd < handle
    %This class groups element of class "element" which describe
    %abcd matrices.
    %By using the absolute distance a complete system is built
    %by adding free space distances
    
    properties
        elements=[];
        elements_pos=[];
    end
    
    methods (Static)
        function ret_q=propagate_gauss(abcd,q)
            % 1,1 - A
            % 1,2 - B
            ret_q=(abcd(1,1).*q+abcd(1,2))./ (abcd(2,1).*q+abcd(2,2));
        end
    end
    
    methods
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
            suptitle('Matrix Optics: Paraxial Beam Propagation');
        end
 
        function [R]=gauss_curvature(obj,z,zr)
            R=z.*(1+(zr./z).^2);  
        end
        
        function [R]=gauss_curvature_q(obj,q)
            R=real(q).*(1+(imag(q)./real(q)).^2);  
        end
        
        function [q,pos, R]=plot_gauss(obj, in_q)

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
            if(nargout==0)
                subplot(3,1,1);
                hold on;
                for i=1:length(obj.elements_pos); xline(obj.elements_pos(i),'--r');end
                plot(pos, real(q));
                xlabel('z (m)');
                ylabel('z(z) (m)');
                subplot(3,1,2);
                stairs(pos, imag(q));
                xlabel('z (m)');
                ylabel('z0(z) (rad)');
                
                subplot(3,1,3);
                plot(pos, R,'*--');
                xlabel('z (m)');
                ylabel('R(z) (m)');
                
                suptitle('Matrix Optics: Gauss Beam Propagation');
            end
        end 
        
    end
end

