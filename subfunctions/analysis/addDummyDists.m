function new_distmat = addDummyDists(distmat,n1)

new_distmat = zeros(20,20);
n2 = size(distmat,1)-n1;

% within target
new_distmat(1:n1,1:n1) = distmat(1:n1,1:n1);
new_distmat(11:10+n2,11:10+n2) = distmat(n1+1:end,n1+1:end);

% between targets
new_distmat(11:10+n2,1:n1) = distmat(n1+1:end,1:n1);
new_distmat(1:n1,11:10+n2) = distmat(1:n1,n1+1:end);

end