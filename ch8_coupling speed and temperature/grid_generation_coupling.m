function [JXYV, JXYP, JMV, JMP, noElement, noNodesV, noNodesP, BE1,BE2,BE3,BE4,BP1,BP2,BP3,BP4] =...
    grid_generation_coupling(H, L, Nx, Ny)
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

%%%%%%% 总单元数和结点数
noElement = Nx * Ny;      %总单元数
noNodesV = (2*Nx+1) * (2*Ny+1);  % 二次单元结点总数 Velocity
noNodesP = (Nx+1) * (Ny+1);      % 线性单元结点总数 Pressure
%%%%%%% 总单元数和结点数

%%%%%%% 结点分布拓扑
chan = zeros(Nx*2+1,Ny*2+1);
chan(1:2:Nx*2+1,1:2:Ny*2+1)= reshape(1:noNodesP,Nx+1,Ny+1);
chan(~chan) = noNodesP+1:noNodesV;
%%%%%%% 结点分布拓扑

%%%%%%% 收敛比例系数
SS = 0.5; %收敛比，定义为区域中最窄部分的高度与最宽部分高度的比值
SL1 = Nx/5*2+1;
SL2 = Nx/5*2;
SL3 = Nx/5*2;
SL4 = Nx/5*2;
SL5 = Nx/5*2;
sl = [ones(1,SL1),(1-(1-SS)/SL2):(SS-1)/SL2:SS,ones(1,SL3)*SS,...
   SS+(1-SS)/SL4:(1-SS)/SL4:1,ones(1,SL5)]';
%%%%%%% 收敛比例系数

%%%%%%% 四边形二次单元JXYV生成
JXYV = zeros(noNodesV, 2);
JXYV(chan(:),1) = reshape(repmat(linspace(0, L, 2*Nx+1),2*Ny+1,1)',noNodesV,1);
JXYV(chan(:),2) = reshape(repmat(linspace(0, H, 2*Ny+1),2*Nx+1,1).*sl,noNodesV,1);
%%%%%%% 四边形二次单元JXYV生成

%%%%%%% 四边形二次单元JMV生成
JMV = zeros(noElement, 9);
tag = 1;
for j = 1:Ny
    for i = 1:Nx
        JMV(tag,:) = reshape(chan(2*i-1:2*i+1,2*j-1:2*j+1), 1, 9);
        tag = tag + 1;
    end
end
%%%%%%% 四边形二次单元JMV生成

%%%%%%% 四边形线性单元JMP和JXYP生成
JMP = [JMV(:,1),JMV(:,3),JMV(:,9),JMV(:,7)];
JXYP = JXYV(1:noNodesP,:);
%%%%%%% 四边形线性单元JMP和JXYP生成

%%%%%%% BP数据生成
[BP1,BP2,BP3,BP4] = deal(chan(:,1), chan(end,:), chan(:,end),chan(1,:));
%%%%%%% BP数据生成

%%%%%%% BE数据生成
BE1 = cat(2, (1:Nx)', ones(Nx,1), zeros(Nx,1), -ones(Nx,1));
BE2 = cat(2, (Nx:Nx:noElement)', 2*ones(Ny,1), ones(Ny,1), zeros(Ny,1));
BE3 = cat(2, (noElement-Nx+1:noElement)', 3*ones(Nx,1), zeros(Nx,1), ones(Nx,1));
BE4 = cat(2, (1:Nx:noElement-Nx+1)', 4*ones(Ny,1), -ones(Ny,1), zeros(Ny,1));
%%%%%%% BE数据生成

%%%%%%% 存储网格数据
% save('mesh_coupling.mat', 'noElement', 'noNodesV', 'noNodesP','JMV','JMP','JXYV',...
%     'JXYP', 'BP1', 'BP2', 'BP3', 'BP4', 'BE1', 'BE2', 'BE3', 'BE4')
%%%%%%% 存储网格数据