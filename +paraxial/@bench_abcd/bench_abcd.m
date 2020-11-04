classdef bench_abcd < handle
    %This class groups element of class "element" which describe
    %abcd matrices.
    %By using the absolute distance a complete system is built
    %by adding free space distances
    
    properties
        elements=[];
        elements_pos=[];
        rl =[0 0]
    end
    
    methods
        function obj = bench_abcd(in_r0, in_alpha0)
            %Initial Beam Parameter
            obj.rl=[in_r0 in_alpha0]';
        end
        
        function add(obj, in_pos,in_element)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.elements=[obj.elements in_element];
            obj.elements_pos    =[obj.elements_pos in_pos];
        end
        
        function plot(obj)
             figure();
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
                    last_rl=obj.rl;
                end

                z= obj.elements_pos(i)-last_z;
                tmp_rl=[1,z;0,1]*last_rl; %free space
                this_rl(:,i)=obj.elements(i).abcd*tmp_rl; 
            end
            
            % Build complete with initial
            rl= [obj.rl, this_rl];
            pos=[0,obj.elements_pos];
            plot(pos, rl(1,:));
            xlabel('z (m)');
            ylabel('r(z) (m)');
            subplot(2,1,2);
            stairs(pos, rl(2,:));
            xlabel('z (m)');
            ylabel('alpha(z) (rad)');
        end
        
    end
end

