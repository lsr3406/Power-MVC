%% 纯文本文件输出
classdef Plain < handle
	properties
		% nothing
	end
	methods

		%% getPowerFlowReport: 返回电力系统潮流计算结果 (TODO string)
		function [result] = getPowerFlowReport(self, steadyState, nodes, generator, branches, config, result)
			
			file_t = fopen(config.documentName,'w');

			fprintf(file_t, '%s\n', '=========================================================================');
			fprintf(file_t, '%s\n', '                          电力系统潮流计算报告');
			fprintf(file_t, '%s\n', '=========================================================================');

			fprintf(file_t, '\n%s%d%s%d%s%d%s%d\n','节点数       PQ / PV / 平衡 / 合计      ',nodes.pqCount,' / ',nodes.pvCount,' / ',nodes.refCount,' / ',length(nodes.id));
			fprintf(file_t, '%s%d%s%d%s%d\n','支路数       普通线路 / 变压器 / 合计      ',branches.lineCount,' / ',branches.transformerCount,' / ',length(branches.id));
			fprintf(file_t, '%s%d%s%d%s%d\n','发电机数       调压厂 / 调频厂 / 合计      ',generator.pvCount,' / ',generator.refCount,' / ',generator.pvCount+generator.refCount);
			fprintf(file_t, '%s%d\n','负荷数                                  ',nodes.loadCount);

			fprintf(file_t, '\n%s%s%d\n', '计算方法: ', config.method);

			if result.status == 101
				fprintf(file_t, '\n%s\n', '迭代次数超过最大值,计算失败');
				fclose(file_t);
				return;
			elseif result.status ~= 1
				fprintf(file_t, '\n%s\n', '未知错误: ', result.status);
				fclose(file_t);
				return;
			end

			fprintf(file_t, '\n%s%d\n', '迭代次数: ',result.it);

			fprintf(file_t, '\n%s\n', '节点电压,功率');
			fprintf(file_t, '\n%s\n', ' id     电压      相角   发电机有功  发电机无功  负荷有功    负荷无功   补偿有功  补偿无功');
			for k = 1:length(steadyState.nodes.id)
				fprintf(file_t, '  %d',steadyState.nodes.id(k));
				fprintf(file_t, '   %6.3f',steadyState.nodes.mag(k),steadyState.nodes.ang(k).*180./pi,steadyState.nodes.Pg(k).*100,steadyState.nodes.Qg(k).*100,steadyState.nodes.Pd(k).*100,steadyState.nodes.Qd(k).*100,steadyState.nodes.Pc(k).*100,steadyState.nodes.Qc(k).*100);
				fprintf(file_t, '\n');
			end
			fprintf(file_t, '%s%6.3f  %6.3f   %6.3f  %6.3f   %6.3f  %6.3f','合计:                     ',sum(steadyState.nodes.Pg).*100,sum(steadyState.nodes.Qg).*100,sum(steadyState.nodes.Pd).*100,sum(steadyState.nodes.Qd).*100,sum(steadyState.nodes.Pc.*100),sum(steadyState.nodes.Qc.*100));
			
			% fprintf(file_t, '\n%s\n', '平衡,PV结点发电功率');
			% fprintf(file_t, '   %s\n','id type  Pg      Qg');

			% for k = sort([steadyState.nodes.Balance_index;steadyState.nodes.PV_index]')

			% 	fprintf(file_t, '   %d   %d  %6.4f  %6.4f\n',steadyState.nodes.id(k),steadyState.nodes.type(k),steadyState.nodes.Pg(k)*100,steadyState.nodes.Qg(k)*100);
			% end
			% fprintf(file_t, ' %s\n','发电及负荷总量');
			% fprintf(file_t, '   %6.4f %6.4f    %6.4f %6.4f',sum(steadyState.nodes.Pg)*100,sum(steadyState.nodes.Qg)*100,sum(steadyState.nodes.Pd)*100,sum(steadyState.nodes.Qd)*100);
			fprintf(file_t, '\n');


			fprintf(file_t, '\n%s\n\n', '线路功率');
			fprintf(file_t, '%s\n', 'id from to    Pij      Qij        Pji      Qji        dP      dQ');
			for k = steadyState.branches.id'
				fprintf(file_t, '  %d',[steadyState.branches.id(k),steadyState.branches.fid(k),steadyState.branches.tid(k)]);
				fprintf(file_t, '   %6.3f',[steadyState.branches.Pij(k),steadyState.branches.Qij(k),steadyState.branches.Pji(k),steadyState.branches.Qji(k),steadyState.branches.dP(k),steadyState.branches.dQ(k)]*100);
				fprintf(file_t, '\n');
			end
			
			fprintf(file_t, '总线损');
			fprintf(file_t, '   %6.3f',sum(steadyState.branches.dP)*100,sum(steadyState.branches.dQ)*100);

			fprintf(file_t, '\n');

			fprintf(file_t, '\n');



			fclose(file_t);

		end

	end
end