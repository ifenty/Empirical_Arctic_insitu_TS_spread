function time_qc = ocean_time_quality_check(time,time_range)
% time_range��time����ʼ��Χ
time_qc = time>=time_range(1,1)&time<time_range(1,2);
end