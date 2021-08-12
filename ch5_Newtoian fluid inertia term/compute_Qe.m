function [Qe1, Qe2] = compute_Qe(JXYe, velocityE, rho, gaussPts, gaussWts);
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


%%%%%%% 初始化Qe1和Qe2
Qe1 = zeros(9, 9);
Qe2 = zeros(9, 9);
%%%%%%% 初始化Qe1和Qe2

xi = gaussPts;
et = gaussPts;

%%%%%%% Qe1和Qe2的数值积分
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
        commonFac = rho * gaussWts(i) * gaussWts(j) * det_Jacobi * (fyV * fyV');
        
        Qe1 = Qe1 + commonFac * velocityE(:,1) * fy_x(1,:);
        Qe2 = Qe2 + commonFac * velocityE(:,2) * fy_x(2,:);
        %%%%%%% Qe1和Qe2单元方程子块计算
    end
end
%%%%%%% Qe1和Qe2的数值积分