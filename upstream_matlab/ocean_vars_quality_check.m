function [lon_qc,lat_qc,time_qc,temp_qc,salt_qc,pres_qc,dens_qc] = ocean_vars_quality_check(lon,lat,time,temp,salt,pres,dens)
% ������γ�ȼ��飬�������ݶȼ��飬������ֵ���飬��ȵߵ��Ժ��ظ��Լ���
% ��ͨ��Ϊ0��ͨ��һ��Ϊ1��ͨ��2��Ϊ2
lon_range = [-180,180];
lat_range = [60,90];
time_range = [datetime(1983,1,1,0,0,0),datetime(2024,1,1,0,0,0)];
depth_std = [0;10;20;30;50;75;100;125;150;200;250;300;400;500;600;700;800;...
    900;1000;1100;1200;1300;1400;1500;1750;2000;2500;3000;3500;4000;4500;5000;...
    5500];
Tmax = [20;20;20;14;14;14;14;14;10;10;10;10;10;10;9;9;9;9;8;8;8;8;8;8;8;8;8;7;7;...
    7;7;7;3];
Tmin = [-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;...
    -2;-2;-2;-2-2;-2;-1.5;-1.5;-1.5;-1.5];
temp_threshold = [Tmax,Tmin,depth_std];
Smax = [40;40;40;40;40;40;38;38;38;38;38;38;37;37;37;37;37;37;37;36;36;36;36;36;...
    36;36;35.5;35.5;35.5;35.5;35.5;35.5;35.5];
Smin = [0;0;0;0;0;0;26;26;26;26;26;30;33;33;33;33;33;33;33;33;33;33;33;33;33;33;...
    33;33;33;33;33;33;33];
salt_threshold = [Smax,Smin,depth_std];
lon_qc = ocean_lon_quality_check(lon,lon_range);
lat_qc = ocean_lat_quality_check(lat,lat_range);
time_qc = ocean_time_quality_check(time,time_range);
temp_qc = ocean_temp_quality_check(temp,pres,temp_threshold);
salt_qc = ocean_salt_quality_check(salt,pres,salt_threshold);
pres_qc = ocean_pres_quality_check(pres);
dens_qc = ocean_dens_quality_check(dens,pres);
end