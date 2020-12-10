function ppp_missed = ppp_miss(ppp,measmodel)
%PPP_MISS performs the misdetection update of a PPP

%the sensor fails to the detect the target
Qd_1 = 1 - measmodel.Pd;

ppp_missed = ppp;
for i = 1:length(ppp_missed)
    %the sensor detects the target, but the target generates 0 measurement
    Qd_2 = measmodel.Pd*(ppp(i).beta/(ppp(i).beta+1))^ppp(i).alpha;
    Qd = Qd_1 + Qd_2;
    ppp_missed(i).w = ppp(i).w + log(Qd);
    ppp_missed(i).beta = 1/(Qd_1/Qd/ppp(i).beta+Qd_2/Qd/(ppp(i).beta+1));
end

end

