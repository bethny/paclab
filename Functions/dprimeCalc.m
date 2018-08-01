function [data,dprime] = dprimeCalc(var,rsp,same,diff)

    hit = 0; miss = 0; falseAlarm = 0; corReject = 0;
    for i = 1:length(var)
        if var(i) && rsp(i)
            hit = hit + 1;
        elseif var(i) && ~rsp(i)
            miss = miss + 1;
        elseif ~var(i) && rsp(i)
            falseAlarm = falseAlarm + 1;
        elseif ~var(i) && ~rsp(i)
            corReject = corReject + 1;
        end
    end
    
    H = hit/diff;
    if H == 0
        H = 0.01;
    end
    F = falseAlarm/same;
    if F == 0
        F = 0.01;
    end
    dprime = norminv(H) - norminv(F);
    
    data.hit = hit;
    data.miss = miss;
    data.falseAlarm = falseAlarm;
    data.corReject = corReject;
end