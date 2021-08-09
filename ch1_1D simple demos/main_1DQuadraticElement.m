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
% d2u/dx2 - 1 = 0, 0<=x<=1
% s.t.     u = 1,   if x = 0
%      du/dx = 2,   if x = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%    网格离散数据     %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noNodes = 11;
noElement = (noNodes-1) / 2;
JX = linspace(0,1,noNodes).';
JM = [1:2:noNodes-1; 2:2:noNodes; 3:2:noNodes].';
JB1 = [1, 1];
JB2 = [noElement, 3, 2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%    网格离散数据     %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 定义数值积分参数
gaussPts = [0.932469514203152,0.661209386466265,0.238619186083197,-0.932469514203152,-0.661209386466265,-0.238619186083197];
gaussWts = [0.171324492379170,0.360761573048139,0.467913934572691,0.171324492379170,0.360761573048139,0.467913934572691];
%%%% 定义数值积分参数

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% stiffness和force的计算 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stiffness = zeros(noNodes);  %初始化K
force = zeros(noNodes,1);  %初始化b
for i=1:noElement 
    JXe = JX(JM(i,:),1);
    stiffness(JM(i,:),JM(i,:)) = stiffness(JM(i,:),JM(i,:)) + compute_stiffnessElement(JXe, gaussPts, gaussWts);
    force(JM(i,:),1) = force(JM(i,:),1) + compute_forceElement(JXe, gaussPts, gaussWts);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% stiffness和force的计算 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Neumann边界条件，代入JB2，计算b和F相加 %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(JB2(:,1)) % 循环次数为JB2的行数
    force(JM(JB2(i,1),JB2(i,2))) = force(JM(JB2(i,1),JB2(i,2))) + JB2(i,3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Neumann边界条件，代入JB2，计算b和F相加 %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% 代入Dirichlet边界条件，求解 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prescribedDof = JB1(:,1);
% free Dof : activeDof
activeDof = setdiff( (1:noNodes)',(prescribedDof) );
force = force - stiffness(:,prescribedDof) * JB1(:,2);

% solution 求解方程 
solution = zeros(noNodes,1);
solution(activeDof) = stiffness(activeDof,activeDof) \ force(activeDof);
solution(prescribedDof) = JB1(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% %%%%% %    代入Dirichlet边界条件，求解     %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  结果对比 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result = [JX,solution]  %结果输出

x = linspace(0,1,101); % 精确解定义域
exactSolution = 1/2*x.^2 + x + 1; % 精确解表达式
plot(JX,solution,'b-*', x, exactSolution, 'r--')  %数值解绘图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% 求解方程，结果对比 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%