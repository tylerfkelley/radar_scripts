function S = flip_phase(x)
m = abs(x);
a = -angle(x);
S = m.*cos(a) + j.*m.*sin(a);
