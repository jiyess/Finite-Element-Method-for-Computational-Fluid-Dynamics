function [JXY,JM,noElement,noNodes,BE1,BE2,BE3,BE4,BP1,BP2,BP3,BP4] = grid_generation_tri(H, L, Nx, Ny)
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
noElement = 2 * Nx * Ny;      % 总单元数
noNodes = (Nx+1) * (Ny+1);  % 结点总数
%%%%%%% 总单元数和结点数

%%%%%%% 结点分布拓扑
AAA = reshape(1:noNodes,Nx+1,Ny+1)';
%%%%%%% 结点分布拓扑

%%%%%%% 三角形单元JXY生成
JXY = zeros(noNodes,2);
JXY(:,1) = repmat(linspace(0,L,Nx+1)',Ny+1,1);
JXY(:,2) = reshape(repmat(linspace(0,H,Ny+1),Nx+1,1),noNodes,1);
%%%%%%% 三角形单元JXY生成

%%%%%%% 三角形单元JM生成
JM = zeros(noElement, 3);
k = 0;
for i=1 : Ny/2
    for j=1 : Nx
        k = k+1;
        JM(k,:) = [AAA(i,j), AAA(i,j+1), AAA(i+1,j)];
        k = k+1;
        JM(k,:) = [AAA(i+1,j+1), AAA(i+1,j), AAA(i,j+1)];
    end
end

for i = Ny/2+1 : Ny
    for j = 1 : Nx
        k = k+1;
        JM(k,:) = [AAA(i+1,j+1), AAA(i,j), AAA(i,j+1)];
        k = k+1;
        JM(k,:) = [AAA(i,j), AAA(i+1,j+1), AAA(i+1,j)];
    end
end
%%%%%%% 三角形单元JM生成

%%%%%%% BP数据生成
BP1 = AAA(1,:);
BP2 = AAA(:,Nx+1);
BP3 = AAA(Ny+1,:);
BP4 = AAA(:,1);
%%%%%%% BP数据生成

%%%%%%% BE数据生成
BE1 = cat(2, (1:2:2*Nx)', ones(Nx,1), zeros(Nx,1), -ones(Nx,1));
BE2 = cat(2, [(2*Nx:2*Nx:Ny*Nx)'; (Ny*Nx+2*Nx-1:2*Nx:2*Ny*Nx)'], 3*ones(Ny,1), ones(Ny,1), zeros(Ny,1));
BE3 = cat(2, (2*Nx*(Ny-1)+2:2:2*Nx*Ny)', 2*ones(Nx,1), zeros(Nx,1), ones(Nx,1));
BE4 = cat(2, [(1:2*Nx:Nx*(Ny-1)+1)';(Nx*Ny+2:2*Nx:2*Nx*Ny)'], 3*ones(Ny,1), -ones(Ny,1), zeros(Ny,1));
%%%%%%% BE数据生成

%%%%%%% 存储网格数据
% save('msh.mat', 'noElement', 'noNodes', 'JM', 'JXY', ...
%     'BP1', 'BP2', 'BP3', 'BP4', 'BE1', 'BE2', 'BE3', 'BE4')
%%%%%%% 存储网格数据