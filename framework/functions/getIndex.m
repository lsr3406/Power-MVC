%% getIndex: 返回第二个数组的元素在第一个数组中的位置
function [index] = getIndex(arr, obj)
	index = [];
	for k = 1:length(obj)
		index = [index, find(arr == obj(k))];
	end
end