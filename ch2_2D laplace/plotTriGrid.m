function plotTriGrid(JM,JXY,numberElement,numberNodes,num_prnt)
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

figure
hold on

triplot(JM, JXY(:,1), JXY(:,2));

if num_prnt == 1
    del_x=0.01;
    del_y=0;%结点号放置的位置的调整
    
    textPosition = [mean(JXY(JM),2)-del_x, mean(JXY(JM+numberNodes),2)-del_y];
    for e = 1:numberElement
        text(textPosition(e,1), textPosition(e,2), int2str(e))
    end
    
    for n = 1:numberNodes
        text(JXY(n,1)-del_x, JXY(n,2)-del_y, ['(',int2str(n),')'])
    end
end

axis equal tight off
title('Triangular Mesh ')
