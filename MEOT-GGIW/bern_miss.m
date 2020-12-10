function [l_missed,bern_missed] = bern_miss(bern,measmodel)
%BERN_MISS performs the misdetection update of a Bernoulli and computes the
%misdetection likelihood

%the sensor fails to the detect the target
Qd_1 = 1 - measmodel.Pd;
%the sensor detects the target, but the target generates 0 measurement
Qd_2 = measmodel.Pd*(bern.beta/(bern.beta+1))^bern.alpha;

%effective misdetection probability
Qd = Qd_1 + Qd_2;

temp = 1-bern.r+bern.r*Qd;
l_missed = log(temp);

bern_missed = bern;
bern_missed.r = bern.r*Qd/temp;
bern_missed.beta = 1/(Qd_1/Qd/bern.beta+Qd_2/Qd/(bern.beta+1));

end

