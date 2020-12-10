function bern_predicted = bern_pred(bern,motionmodel)
%BERN_PRED performs the prediction of a Bernoulli component

bern_predicted.r = bern.r*motionmodel.Ps;

bern_predicted.xr = motionmodel.Ar*bern.xr;
bern_predicted.Cr = motionmodel.Ar*bern.Cr*motionmodel.Ar'+motionmodel.Cwr;

bern_predicted.v = 6 + exp(-motionmodel.Ts/motionmodel.tau)*(bern.v-6);
bern_predicted.V = exp(-motionmodel.Ts/motionmodel.tau)*bern.V;

bern_predicted.alpha = bern.alpha/motionmodel.eta;
bern_predicted.beta = bern.beta/motionmodel.eta;

end

