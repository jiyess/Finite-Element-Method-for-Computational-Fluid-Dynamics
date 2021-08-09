function stiffnessElement = compute_stiffnessElement(JXe, gaussPts, gaussWts)

%%% 初始化 stiffnessElement
stiffnessElement = zeros(length(JXe));

for i = 1 : length(gaussWts)
    %%%% 插值函数对kesi的倒数
    dfy_dxi=[gaussPts(i)-1/2;  -2*gaussPts(i); gaussPts(i)+1/2 ];      % 式（1.47）
    %%%% 插值函数对kesi的倒数
    
    %%%%  Jacobi相关
    dx_dxi = dfy_dxi' * JXe;  % 式（1.57）
    J = dx_dxi;               %式（1.55）
    dfy_dx = dfy_dxi / J;   % 式（1.54）
    %%%%  Jacobi相关
    
    %%%%  计算stiffnessElement
    stiffnessElement = stiffnessElement + gaussWts(i) * (dfy_dx * dfy_dx') * det(J);   % 式（1.62）
    %%%%  计算stiffnessElement
    
end