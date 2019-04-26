% encoding: utf-8
% @author: siru
% @create on: 2019-03-05 15:52:03
% @update on: 2019-04-26 14:44:15

% 数据类型检验器
classdef Validate < handle

properties (Constant)

    % pf
    pfDefaultConfig = {{'method', 'NR'}, {'n_iters_max', 50}, {'epsilon', '1e-6'}, {'start', 'flat'}};
    pfMethods = {'NR', 'FD', 'FDBX', 'FDXB'};
    pfStart = {'flat', 'plain', 'dc'};

    % opf
    opfDefaultConfig = {{'n_iters_max', 50}, {'epsilon', 1e-6}, {'sigma', 1e-2}};

    % ld
    ldDefaultConfig = {{'dae', 'euler'}, {'net', 'default'}, {'dt', 0.001}, {'time', 2}};
    ldDae = {'euler'};
end

methods (Static)


end

end
