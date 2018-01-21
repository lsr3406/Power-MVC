%% gamma2pi: 将变压器的γ形等效电路参数转化成π形等效电路参数
function [from_o, to_o, Z_o, Y1_o, Y2_o] = gamma2pi(nodes,from,to,Z,Y,K)
	
	from_o = from;
	to_o = to;

	Z_o = K.*Z;
	Y1_o = (1-K)./Z./K.^2 + Y;
	Y2_o = (K-1)./Z./K;

end