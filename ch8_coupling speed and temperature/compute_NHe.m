function NHe = compute_NHe(JXYe, velocityE, muE, gaussPts, gaussWts)
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

%%%%%%% 初始化NHe
NHe=zeros(9,1);
%%%%%%% 初始化NHe

xi = gaussPts;
et = gaussPts;

for i = 1:length(gaussWts)
    for j = 1:length(gaussWts)
        
        %%%%%%%% 速度插值函数对xi和et的导数
        fyV = [1/4*xi(i)*(-1)*(xi(i)-1)*((-1)-1);
            1/2*(-1)*(1-xi(i)^2) *((-1)-1);
            1/4*xi(i)*(-1)*(xi(i)+1)*((-1)-1);
            1/2*xi(i)*( xi(i)-1)*(1-(-1)^2);
            (1-xi(i)^2)*(1-(-1)^2);
            1/2*xi(i)*( xi(i)+1)*(1-(-1)^2);
            1/4*xi(i)*(-1)*(xi(i)-1)*((-1)+1);
            1/2*(-1)*(1-xi(i)^2) *((-1)+1);
            1/4*xi(i)*(-1)*(xi(i)+1)*((-1)+1)];
        
        fyV_xi = [1/4*et(j)*(xi(i)-1)*(et(j)-1) + 1/4*xi(i)*et(j)*(et(j)-1);
            -et(j)*xi(i)*(et(j)-1);
            1/4*et(j)*(xi(i)+1)*(et(j)-1)+1/4*xi(i)*et(j)*(et(j)-1);
            1/2*(xi(i)-1)*(1-et(j)^2)+1/2*xi(i)*(1-et(j)^2);
            -2*xi(i)*(1-et(j)^2);
            1/2*(xi(i)+1)*(1-et(j)^2)+1/2*xi(i)*(1-et(j)^2);
            1/4*et(j)*(xi(i)-1)*(et(j)+1)+1/4*xi(i)*et(j)*(et(j)+1);
            -et(j)*xi(i)*(et(j)+1);
            1/4*et(j)*(xi(i)+1)*(et(j)+1)+1/4*xi(i)*et(j)*(et(j)+1)];
        
        fyV_et = [1/4*xi(i)*(xi(i)-1)*(et(j)-1)+1/4*xi(i)*et(j)*(xi(i)-1);
            1/2*(1-xi(i)^2)*(et(j)-1)+1/2*et(j)*(1-xi(i)^2);
            1/4*xi(i)*(xi(i)+1)*(et(j)-1)+1/4*xi(i)*et(j)*(xi(i)+1);
            -xi(i)*et(j)*(xi(i)-1);
            -2*(1-xi(i)^2)*et(j);
            -xi(i)*et(j)*(xi(i)+1);
            1/4*xi(i)*(xi(i)-1)*(et(j)+1)+1/4*xi(i)*et(j)*(xi(i)-1);
            1/2*(1-xi(i)^2)*(et(j)+1)+1/2*et(j)*(1-xi(i)^2);
            1/4*xi(i)*(xi(i)+1)*(et(j)+1)+1/4*xi(i)*et(j)*(xi(i)+1)];
        %%%%%%%% 速度插值函数对xi和et的导数
        
        %%%%%%%% Jacobi相关计算
        Jacobi = [fyV_xi, fyV_et]' * JXYe;
        fy_x = Jacobi \ [fyV_xi, fyV_et]';
        
        det_Jacobi = det(Jacobi);
        %%%%%%%% Jacobi相关计算

       %%%%%%%% 黏度计算
       mu = fyV' * muE;
       %%%%%%%% 黏度计算

       %%%%%%% NHe单元方程子块计算
       NHe = NHe + gaussWts(i)*gaussWts(j)*mu*fyV*((2*fy_x(1,:)*velocityE(:,1)*fy_x(1,:)*velocityE(:,1))+...
           (2*fy_x(2,:)*velocityE(:,2)*fy_x(2,:)*velocityE(:,2))+(fy_x(2,:)*velocityE(:,1)+fy_x(1,:)*velocityE(:,2))^2)*det_Jacobi;
       %%%%%%% NHe单元方程子块计算
    end
end