function [Fe1, Fe2] = compute_forceElement(JXYe,JBPe,gaussPts,gaussWts)
%    see Bi Chao, "Finite Element Method for Computational Fluid
%    Dynamics and Its Detailed Programming"
%    中文：毕超《计算流体力学有限元方法及其编程详解》
%    
%    Copyright (C) 2021 Ye Ji, Dalian University of Technology
%    
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%    
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%    
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

%%%%%%% 初始化Fe1和Fe2
Fe1=zeros(9,1);
Fe2=zeros(9,1);
%%%%%%% 初始化Fe1和Fe2

xi = gaussPts;
et = gaussPts;

switch JBPe(1,2)
    case 1
        for i = 1:length(gaussWts)
            
            fyV = shapeFunction(xi(i), -1 , 'quadratic');
            fyP = shapeFunction(xi(i), -1 , 'linear');
            
            P=[JBPe(1,5);JBPe(1,6);0;0];
            
            le = norm(JXYe(1,:)-JXYe(3,:),2);
            
            commonFac = 0.5 * gw(i) * fyP' * P * le * fyV;
            Fe1 = Fe1 + commonFac * JBPe(1,3);
            Fe2 = Fe2 + commonFac * JBPe(1,4);
        end
        
    case 2
        
        for j = 1:length(gaussWts)
            
            fyV = shapeFunction(1, et(j) , 'quadratic');
            fyP = shapeFunction(1, et(j) , 'linear');
            
            P=[0;JBPe(1,5);JBPe(1,6);0];
            
            le = norm(JXYe(3,:)-JXYe(9,:),2);
            
            commonFac = 0.5 * gaussWts(j) * fyP' * P * le * fyV;
            Fe1 = Fe1 + commonFac * JBPe(1,3);
            Fe2 = Fe2 + commonFac * JBPe(1,4);
        end
        
    case 3
        
        for i = 1:length(gaussWts)
            
            fyV = shapeFunction(xi(i), 1 , 'quadratic');
            fyP = shapeFunction(xi(i), 1 , 'linear');
            
            P=[0;0;JBPe(1,5);JBPe(1,6)];
            
            le = norm(JXYe(7,:)-JXYe(9,:),2);
            
            commonFac = 0.5 * gaussWts(i) * fyP' * P * le * fyV;
            Fe1 = Fe1 + commonFac * JBPe(1,3);
            Fe2 = Fe2 + commonFac * JBPe(1,4);
        end
        
    case 4
        for j = 1:length(gaussWts)
            
            fyV = shapeFunction(-1, et(j) , 'quadratic');
            fyP = shapeFunction(-1, et(j) , 'linear');
            
            P=[JBPe(1,6);0;0;JBPe(1,5)];
            
            le = norm(JXYe(7,:)-JXYe(1,:),2);
            
            commonFac = 0.5 * gaussWts(j) * fyP' * P * le * fyV;
            Fe1 = Fe1 + commonFac * JBPe(1,3);
            Fe2 = Fe2 + commonFac * JBPe(1,4);
        end
        
end