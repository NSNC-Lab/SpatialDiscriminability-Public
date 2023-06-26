function [ll,ul,av] = sem(mat)

mat(isnan(mat)) = [];
av = mean(mat);

if isvector(mat)
    ll=av-(std(mat)/sqrt(length(mat)));
    ul=av+(std(mat)/sqrt(length(mat)));
else
    ll=av-(std(mat)/sqrt(size(mat,1)));
    ul=av+(std(mat)/sqrt(size(mat,1)));
end

end
