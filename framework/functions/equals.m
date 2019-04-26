% encoding: utf-8
% @author: siru
% @create on: 2019-04-26 14:31:23
% @update on: 2019-04-26 14:46:40
%% equals: 判断两个矩阵在数值上是否满足给定的要求
function [res] = equals(mat1, mat2, tolerance)
    if nargin == 2
        tolerance = 1e-8;
    end
    res = all(all( abs(mat1 - mat2) <= tolerance ));
end
