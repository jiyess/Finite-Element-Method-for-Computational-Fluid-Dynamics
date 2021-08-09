function forceElement = compute_forceElement(JXe, gaussPts, gaussWts)

%%% 初始化forceElement
forceElement = zeros(length(JXe),1);
%%% 初始化forceElement

for i = 1 : length(gaussWts)
%%%% 插值函数及其倒数
   fy=[1/2*gaussPts(i)*(gaussPts(i)-1); (1-gaussPts(i))*(1+gaussPts(i)); 1/2*gaussPts(i)*(gaussPts(i)+1)];   % 式（1.46）
   dfy_dxi=[gaussPts(i)-1/2; -2*gaussPts(i); gaussPts(i)+1/2 ]; % 式（1.47）
%%%% 插值函数及其倒数

%%%%  Jacobi相关
    dx_dxi = dfy_dxi'*JXe;   % 式（1.57）
    J = dx_dxi;    %   式（1.55）
%%%%  Jacobi相关

%%%%  计算forceElement
    forceElement = forceElement - gaussWts(i) * fy * det(J);  % 式（1.63）
%%%%  计算forceElement

end