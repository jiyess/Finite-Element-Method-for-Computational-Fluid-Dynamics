function [JXYV, JXYP, JMV, JMP, noElement, noNodesV, noNodesP, BE1,BE2,BE3,BE4,BP1,BP2,BP3,BP4] = grid_generation_inlet_outlet(H,L,Nx,Ny)
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

%%%%%%%  四边形二次单元JXYV生成
JXYV = zeros(noNodesV, 2);
JXYV(chan(:),1) = reshape(repmat(linspace(0, L, 2*Nx+1),2*Ny+1,1)',noNodesV,1);
JXYV(chan(:),2) = reshape(repmat(linspace(0, H, 2*Ny+1),2*Nx+1,1),noNodesV,1);
%%%%%%%  四边形二次单元JXYV生成

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

% 这里开始与第三章网格生成代码不同
%%%%%%% BP数据生成
[BP1,BP2,BP3,BP4] = deal([chan(:,1);chan(end,1:2*Ny/4*3+1)'],...
    chan(end,2*Ny/4*3+1:end), [chan(:,end); chan(1,2*Ny/4*1+1:end)'],...
    chan(1,1:2*Ny/4*1+1));
%%%%%%% BP数据生成

%%%%%%% BE数据生成
BE_down = [(1:Nx)', ones(Nx,1), zeros(Nx,1), -ones(Nx,1)];
BE_right = [(Nx:Nx:noElement)', 2*ones(Ny,1), ones(Ny,1), zeros(Ny,1)];
BE_up = [(noElement-Nx+1:noElement)', 3*ones(Nx,1), zeros(Nx,1), ones(Nx,1)];
BE_left = [(1:Nx:noElement-Nx+1)', 4*ones(Ny,1), -ones(Ny,1), zeros(Ny,1)];

BE1 = cat(1, BE_down, BE_right(1:Ny/4*3,:));
BE2 = BE_right(Ny/4*3+1:Ny,:);
BE3 = cat(1, BE_up, BE_left(Ny/4+1:Ny,:));
BE4 = BE_left(1:Ny/4,:);
%%%%%%% BE数据生成

%%%%%%% 存储网格数据
% save('mesh_inlet_outlet.mat', 'noElement', 'noNodesV', 'noNodesP','JMV','JMP','JXYV',...
%     'JXYP', 'BP1', 'BP2', 'BP3', 'BP4', 'BE1', 'BE2', 'BE3', 'BE4')
%%%%%%% 存储网格数据
