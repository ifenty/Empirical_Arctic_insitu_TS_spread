%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% depth_std = [0;10;20;30;50;75;100;125;150;200;250;300;400;500;600;700;800;...
%     900;1000;1100;1200;1300;1400;1500;1750;2000;2500;3000;3500;4000;4500;5000;...
%     5500];
% Tmax = [20;20;20;14;14;14;14;14;10;10;10;10;10;10;9;9;9;9;8;8;8;8;8;8;8;8;8;7;7;...
%     7;7;7;3];
% Tmin = [-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;...
%     -2;-2;-2;-2-2;-2;-1.5;-1.5;-1.5;-1.5];
% temp_threshold = [Tmax,Tmin,depth_std];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function temp_qc = ocean_temp_quality_check(temp,pres,temp_threshold)
n = size(temp);
len = n(1,2);
wid = n(1,1);
temp_qc = zeros(wid,len);
temp_qc(isnan(temp))=nan;
Tmax = temp_threshold(:,1);
Tmin = temp_threshold(:,2);
depth_std = temp_threshold(:,3);
parfor i = 1:len
    p = pres(:,i);
    t = temp(:,i);% 每个剖面
    local_qc = zeros(wid, 1);
    for j = 1:wid% 每个剖面的每个点
        delta = abs(p(j)-depth_std);
        idx = find(delta == min(delta));
        
        if t(j)>=Tmin(idx)&t(j)<=Tmax(idx)
            local_qc(j) = local_qc(j)+1;
        end
    end
    if sum(~isnan(t))~=1
        for k = 1:wid-1%每个剖面的每个点
            for m = k+1:wid
                diff_value = abs(p(m) - p(k));  % 计算差值
                if diff_value >= 3&local_qc(k)==1  % 如果差值大于等于3，记录该差值并跳出循环
                    Gdt = (t(m)-t(k))./diff_value;
                    if Gdt>=0&&Gdt<=0.3
                        local_qc(k) = local_qc(k)+1;
                    elseif Gdt<0&&abs(Gdt)<=0.7
                        local_qc(k) = local_qc(k)+1;
                    end
                    break
                end
            end
        end
        ind0 = max(find(local_qc~=0,1,'last')-1,1);
        if  local_qc(ind0)==2
            local_qc(ind0+1) = local_qc(ind0+1)+1;
        end
%         if  local_qc(k)==2
%             local_qc(k+1) = local_qc(k+1)+1;
%         end
    else
        ind = find(~isnan(t));
        if local_qc(ind) == 1
            local_qc(ind) = local_qc(ind) + 1;
        end
    end
    temp_qc(:,i) = local_qc;
end