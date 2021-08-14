function CDBe = compute_CDBe(JXYe, JBT2e,  gaussPts, gaussWts)
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

%%%%%%% 初始化CDBe
CDBe=zeros(9,1);
%%%%%%% 初始化CDBe

xi = gaussPts;
et = gaussPts;

%%%%%%% 提取JBT2数据
switch JBT2e(1,2)
    
    case 1
        for i = 1:length(gaussWts)
            %%%%%%%% 速度插值函数对xi和et的导数
            fy = shapeFunction(xi(i), -1 , 'quadratic');
            %%%%%%%% 速度插值函数对xi和et的导数
            
            qb = [JBT2e(1,5:7), 0, 0, 0, 0, 0, 0]';
            q  = fy' * qb;
            
            lb = norm(JXYe(1,:)-JXYe(3,:));
            CDBe = CDBe + gaussWts(i) * q * fy * lb;
        end
        
    case 2
        for j = 1:length(gaussWts)
            
            %%%%%%%% 速度插值函数对xi和et的导数
            fy = shapeFunction(1, et(j) , 'quadratic');
            %%%%%%%% 速度插值函数对xi和et的导数
            
            qb=[0, 0, JBT2e(1,5), 0, 0, JBT2e(1,6), 0, 0, JBT2e(1,7)]';
            q=fy'*qb;
            
            lb = norm(JXYe(3,:)-JXYe(9,:));
            CDBe = CDBe + gaussWts(j) * q * fy * lb;
        end
        
    case 3
        
        for i = 1:length(gaussWts)
            %%%%%%%% 速度插值函数对xi和et的导数
            fy = shapeFunction(xi(i), 1 , 'quadratic');
            %%%%%%%% 速度插值函数对xi和et的导数
            
            qb = [0, 0, 0, 0, 0, 0, JBT2e(1,7:-1:5)]';
            q = fy' * qb;
            
            lb = norm(JXYe(7,:)-JXYe(9,:));
            CDBe = CDBe + gaussWts(i) * q * fy * lb;
        end
        
    case 4
        for j = 1:length(gaussWts)
            
            %%%%%%%% 速度插值函数对xi和et的导数
            fy = shapeFunction(-1, et(j) , 'quadratic');
            %%%%%%%% 速度插值函数对xi和et的导数
            
            qb=[JBT2e(1,7), 0, 0, JBT2e(1,6), 0, 0, JBT2e(1,5), 0, 0]';
            q=fy'*qb;
            
            lb = norm(JXYe(1,:)-JXYe(7,:));
            CDBe = CDBe + gaussWts(j)*q*fy*lb;
        end
end
