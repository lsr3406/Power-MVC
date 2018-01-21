%% gamma2pi: 将变压器的 γ 形等效电路参数转化成 π 形等效电路参数, 按照 matpower 手册的定义
% @param Z: 变压器原始参数-短路阻抗标么值
% @param Y: 变压器原始参数-励磁导纳标么值
% @param K: 变压器非标准变比, 
% @return Z: π 型等效电路线路阻抗(阻抗)
% @return Y1: π 型等效电路一次侧导纳(导纳)
% @return Y2: π 型等效电路二次侧导纳(导纳)
function [Z, Y1, Y2] = gamma2pi(Z,Y,K)
	Y1 = (1-K)./Z./K.^2 + Y./2;
	Y2 = (K-1)./Z./K + Y./2;
	Z = K.*Z;
end
