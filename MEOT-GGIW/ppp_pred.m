function ppp_predicted = ppp_pred(ppp,motionmodel,birthmodel)
%PPP_PRED performs the prediction of a PPP

ppp_predicted = ppp;

for i = 1:length(ppp)
    ppp_predicted(i).w = ppp(i).w + log(motionmodel.Ps);
    
    ppp_predicted(i).xr = motionmodel.Ar*ppp(i).xr;
    ppp_predicted(i).Cr = motionmodel.Ar*ppp(i).Cr*motionmodel.Ar'+motionmodel.Cwr;
    
    ppp_predicted(i).v = 6 + exp(-motionmodel.Ts/motionmodel.tau)*(ppp(i).v-6);
    ppp_predicted(i).V = exp(-motionmodel.Ts/motionmodel.tau)*ppp(i).V;

    ppp_predicted(i).alpha = ppp(i).alpha/motionmodel.eta;
    ppp_predicted(i).beta = ppp(i).beta/motionmodel.eta;
end

ppp_predicted = [ppp_predicted;birthmodel];

end

