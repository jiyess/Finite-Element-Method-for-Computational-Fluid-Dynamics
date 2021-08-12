function [Ge1,Ge2] = compute_Ge(JXYe, rho, gaussPts, gaussWts)
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

%%%%%%% 初始化Ge1和Ge2
Ge1 = zeros(9,1);
Ge2 = zeros(9,1);
%%%%%%% 初始化Ge1和Ge2

xi = gaussPts;
et = gaussPts;

%%%%%%% Ge1和Ge2的数值积分
for i = 1 : length(gaussWts)
    for j = 1 : length(gaussWts)
        
        %%%%%%%% 速度插值函数对xi和et的导数
        [fyV, fyV_xi, fyV_et] = shapeFunction(xi(i), et(j) , 'quadratic');
        %%%%%%%% 速度插值函数对xi和et的导数
        
        %%%%%%%% Jacobi相关计算
        Jacobi = [fyV_xi, fyV_et]' * JXYe;
        fy_x = Jacobi \ [fyV_xi, fyV_et]';
        
        det_Jacobi = det(Jacobi);
        %%%%%%%% Jacobi相关计算
        
        commonFac = gaussWts(i) * gaussWts(j) * det_Jacobi;
        
        %%%%%%% Ge1和Ge2单元方程子块计算
        Ge1 = Ge1 + rho * commonFac * fyV * fy_x(1,:) * fyV;
        Ge2 = Ge2 + rho * commonFac * fyV * fy_x(2,:) * fyV;
        %%%%%%% Ge1和Ge2单元方程子块计算
    end
end
%%%%%%% Ge1和Ge2的数值积分