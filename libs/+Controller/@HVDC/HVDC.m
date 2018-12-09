% encoding: utf-8
% @author: siru
% @create on: 2018-10-05 09:56:42
% @update on: 2018-10-20 17:47:07

classdef HVDC < handle

properties
	% nothing
end

methods

	%% 测试 1
	function sensitivityAnalysis(self)

		% 交流线电压
		data.Er = 175000;
		data.Ei = 160000;

		% 整流侧变压器变比
		data.TrList = 1.05:-0.001:0.95;
		data.TrList = [data.TrList(1:end-1), fliplr(data.TrList)];

		% 实验结果
		res.operatingMode = [];
		res.alpha = [];
		res.beta = [];
		res.gamma = [];
		res.Vdr0 = [];
		res.Vdi0 = [];
		res.Vdr = [];
		res.Vdi = [];
		res.Id = [];
		res.Pr = [];
		res.Pi = [];
		res.Qr = [];
		res.Qi = [];
		res.phir = [];
		res.phii = [];

		% 开始实验, 创建模型, 获取数据
		hvdc = Model.HVDC();
		hvdc.init(dcm_test());

		for k = 1:length(data.TrList)

		    hvdc.Tr = data.TrList(k);  % 调整变比参数
		    hvdc.render(data.Er, data.Ei, 1e-3);  % 测试

		    % 保存数据
		    res.operatingMode = [res.operatingMode, hvdc.operatingMode];
		    res.alpha = [res.alpha, hvdc.alpha];
		    res.beta = [res.beta, hvdc.beta];
		    res.gamma = [res.gamma, hvdc.gamma];
		    res.Vdr0 = [res.Vdr0, hvdc.Vdr0];
		    res.Vdi0 = [res.Vdi0, hvdc.Vdi0];
		    res.Vdr = [res.Vdr, hvdc.Vdr];
		    res.Vdi = [res.Vdi, hvdc.Vdi];
		    res.Id = [res.Id, hvdc.Id];
		    res.Pr = [res.Pr, hvdc.Pr];
		    res.Pi = [res.Pi, hvdc.Pi];
		    res.Qr = [res.Qr, hvdc.Qr];
		    res.Qi = [res.Qi, hvdc.Qi];
		    res.phir = [res.phir, hvdc.phir];
		    res.phii = [res.phii, hvdc.phii];

		end

		figure(1);
		hold on;
		grid on;
		plot(data.TrList, res.alpha.*180./pi, 'LineWidth', 1.5, 'Color', [244, 67, 54]./256);
		plot(data.TrList, res.beta.*180./pi, 'LineWidth', 1.5, 'Color', [76, 175, 80]./256);
		plot(data.TrList, res.gamma.*180./pi, 'LineWidth', 1.5, 'Color', [33, 150, 243]./256);
		set(get(gca,'XLabel'),'String','整流站变压器变比 T_r (pu)');
		set(get(gca,'YLabel'),'String','相应控制角 (deg)');
		title('直流系统控制角随整流站交流电压的变化');
		legend('整流站 \alpha','逆变站 \beta','逆变站 \gamma');

		figure(2);
		hold on;
		grid on;
		plot(data.TrList, res.Vdr0 .* 1e-3, 'LineWidth', 1.5, 'Color', [244, 67, 54]./512);
		plot(data.TrList, res.Vdi0 .* 1e-3, 'LineWidth', 1.5, 'Color', [76, 175, 80]./512);
		plot(data.TrList, res.Vdr .* 1e-3, 'LineWidth', 1.5, 'Color', [244, 67, 54]./256);
		plot(data.TrList, res.Vdi .* 1e-3, 'LineWidth', 1.5, 'Color', [76, 175, 80]./256);
		set(get(gca,'XLabel'),'String','整流站变压器变比 T_r (pu)');
		set(get(gca,'YLabel'),'String','直流电压 (kV)');
		title('直流系统电压随整流站交流电压的变化');
		legend('整流站 V_{dr0}','逆变站 V_{di0}', '整流站 V_{dr}','逆变站 V_{di}');

		figure(3);
		hold on;
		grid on;
		plot(data.TrList, res.operatingMode, 'LineWidth', 1.5, 'Color', [244, 67, 54]./256);
		set(get(gca,'XLabel'),'String','整流站变压器变比 T_r (pu)');
		set(get(gca,'YLabel'),'String','运行方式');
		title('直流系统运行方式随整流站交流电压的变化');

		figure(4);
		hold on;
		grid on;
		plot(data.TrList, res.Id, 'LineWidth', 1.5, 'Color', [244, 67, 54]./256);
		set(get(gca,'XLabel'),'String','整流站变压器变比 T_r (pu)');
		set(get(gca,'YLabel'),'String','直流电流 (A)');
		title('直流系统电流随整流站交流电压的变化');

		figure(5);
		hold on;
		grid on;
		plot(data.TrList, res.Pr .* 1e-6, 'LineWidth', 1.5, 'Color', [244, 67, 54]./512);
		plot(data.TrList, res.Pi .* 1e-6, 'LineWidth', 1.5, 'Color', [76, 175, 80]./512);
		plot(data.TrList, res.Qr .* 1e-6, 'LineWidth', 1.5, 'Color', [244, 67, 54]./256);
		plot(data.TrList, res.Qi .* 1e-6, 'LineWidth', 1.5, 'Color', [76, 175, 80]./256);
		set(get(gca,'XLabel'),'String','整流站变压器变比 T_r (pu)');
		set(get(gca,'YLabel'),'String','功率 (MW / MVar)');
		title('直流系统功率随整流站交流电压的变化');
		legend('整流站 P_{r}','逆变站 P_{i}', '整流站 Q_{r}','逆变站 Q_{i}');

		figure(6);
		hold on;
		grid on;
		plot(data.TrList, cos(res.phir), 'LineWidth', 1.5, 'Color', [244, 67, 54]./256);
		plot(data.TrList, cos(res.phii), 'LineWidth', 1.5, 'Color', [76, 175, 80]./256);
		set(get(gca,'XLabel'),'String','整流站变压器变比 T_r (pu)');
		set(get(gca,'YLabel'),'String','功率因数');
		title('直流系统功率因数随整流站交流电压的变化');
		legend('整流站 cos\varphi_{r}','逆变站 cos\varphi_{i}');

	end
	
end

end

