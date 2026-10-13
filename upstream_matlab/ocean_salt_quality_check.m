function salt_qc = ocean_salt_quality_check(salt,pres,salt_threshold)
n = size(salt);
len = n(1,2);
wid = n(1,1);
salt_qc = zeros(wid,len);
salt_qc(isnan(salt))=nan;
Smax = salt_threshold(:,1);
Smin = salt_threshold(:,2);
depth_std = salt_threshold(:,3);
parfor i = 1:len
    p = pres(:,i);
    s = salt(:,i);
    local_qc = zeros(wid, 1);
    for j = 1:wid
        delta = abs(p(j)-depth_std);
        idx = find(delta == min(delta));
        if s(j)>=Smin(idx)&s(j)<=Smax(idx)
            local_qc(j) = local_qc(j)+1;
        end
    end
    if sum(~isnan(s))~=1
        for k = 1:wid-1
            for m = k+1:wid
                diff_value = abs(p(m) - p(k));  % 计算差值
                if diff_value >= 3&local_qc(k)==1   % 如果差值大于等于3，记录该差值并跳出循环
                    Gds = (s(m)-s(k))./diff_value;
                    if Gds>0&&Gds<9
                        local_qc(k) = local_qc(k)+1;
                    elseif Gds<0&&abs(Gds)<=0.05
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
    else
        ind = find(~isnan(s));
        if local_qc(ind) == 1
            local_qc(ind) = local_qc(ind) + 1;
        end
    end
    salt_qc(:,i) = local_qc;
end