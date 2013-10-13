function [w, h] = nmf_mu(a, k, maxiter)

[m, n] = size(a);

w = rand(m, k);
h = rand(k, n);

for i = 1:maxiter
	h = h.*(w'*a)./(w'*w*h + 1e-9);
	w = w.*(a*h')./(w*h*h' + 1e-9);
end
 



