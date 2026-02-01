function t1k = Wilcoxon(allSamples, testsamples)
t1k=0;
N=length(allSamples);
n=length(testsamples);
m=N-n;

for i=1:length(allSamples)
    if isempty(find(testsamples==allSamples(i), 1))==0
        t1k = t1k + i;
    end
end
mu=n*(N+1)/2;
sigma=m*n*(N+1)/12;
t1k=(t1k-mu)^2/sigma;
end