%    see Bi Chao, "Finite Element Method for Computational Fluid
%    Dynamics and Its Detailed Programming"
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

clc;clear;close all
% du/dx - 1 = 0, 0 <= x <= 1
% s.t. u = 0, if x = 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  网格离散数据 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noNodes = 6;
noElement = noNodes-1;
JX = linspace(0,1,noNodes).';
JM = [1:noNodes-1; 2:noNodes].';
JB1 = [1,0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  网格离散数据 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  计算单元刚度阵 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stiffnessElement = [-1/2 1/2;   % 见式（1.23）
    -1/2 1/2];   %本例中由于采用等距网格划分所有单元ke和be一致
delta_x = 1/noElement;
forceElement = delta_x * [1/2;1/2];  % 见式（1.23）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  计算单元刚度阵 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  单刚组装总刚 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stiffness = zeros(noNodes);  % 初始化K，结点总数为N
force = zeros(noNodes, 1);   % 初始化b
for i = 1:noElement
    stiffness(JM(i,:),JM(i,:)) = stiffness(JM(i,:),JM(i,:)) + stiffnessElement;
    force(JM(i,:),1) = force(JM(i,:),1) + forceElement;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  单刚组装总刚 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% 代入Dirichlet边界条件，求解 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prescribedDof = JB1(:,1);
% free Dof : activeDof
activeDof = setdiff( (1:noNodes)',(prescribedDof));
force = force - stiffness(:,prescribedDof) * JB1(:,2);

% solution 求解方程 
solution = zeros(noNodes,1);
solution(activeDof) = stiffness(activeDof,activeDof)\force(activeDof);
solution(prescribedDof) = JB1(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% 代入Dirichlet边界条件，求解 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = [JX, solution] % 显示结果

figure
plot(JX, solution, 'b-*', JX, JX,'r--')
legend('FEM solution','exact solution')
