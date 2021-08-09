function P9 = P4toP9(p4, JMV)
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

P9 = zeros(9,1);
for i=1:length(JMV(:,1))
    %%%%%%% 构建单元压力
    pp=[p4(JMV(i,1)),p4(JMV(i,3)),p4(JMV(i,9)),p4(JMV(i,7))];
    %%%%%%% 构建单元压力
    
    %%%%%%% 结点1压力
    P9(JMV(i,1))=p4(JMV(i,1));
    %%%%%%% 结点1压力
    
    %%%%%%% 结点2压力
    fyP = shapeFunction(0, -1 , 'linear');
    P9(JMV(i,2))=pp*fyP;
    %%%%%%% 结点2压力
    
    %%%%%%% 结点3压力
    P9(JMV(i,3))=p4(JMV(i,3));
    %%%%%%% 结点3压力
    
    %%%%%%% 结点4压力
    fyP = shapeFunction(-1, 0, 'linear');
    P9(JMV(i,4)) = pp * fyP;
    %%%%%%% 结点4压力
    
    %%%%%%% 结点5压力
    fyP = shapeFunction(0, 0, 'linear');
    P9(JMV(i,5)) = pp * fyP;
    %%%%%%% 结点5压力
    
    %%%%%%% 结点6压力
    fyP = shapeFunction(1, 0, 'linear');
    P9(JMV(i,6)) = pp * fyP;
    %%%%%%% 结点6压力
    
    %%%%%%% 结点7压力
    P9(JMV(i,7)) = p4(JMV(i,7));
    %%%%%%% 结点7压力
    
    %%%%%%% 结点8压力
    fyP = shapeFunction(0, 1, 'linear');
    P9(JMV(i,8))=pp*fyP;
    %%%%%%% 结点8压力
    
    %%%%%%% 结点9压力
    P9(JMV(i,9))=p4(JMV(i,9));
    %%%%%%% 结点9压力
end