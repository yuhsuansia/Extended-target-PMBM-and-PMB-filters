function [l_upd,bern_updated] = bern_upd(bern,z,measmodel)
%BERN_UPD performs the measurement update of a Bernoulli and computes the
%measurement likelihood

[l_upd,bern_updated] = ggiw(bern,z,measmodel);
bern_updated.r = 1;
l_upd = l_upd + log(bern.r);

end

