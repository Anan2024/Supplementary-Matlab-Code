function t2k = AB(allSamples, testsamples)
t2k=0;
N=length(allSamples);
n=length(testsamples);
m=N-n;
for i=1:length(allSamples)
    if isempty(find(testsamples==allSamples(i), 1))==0
        t2k = t2k + abs(i-(N+1)/2);
    end
end
if mod(N,2)==0;
    mu=n*N/4;
    sigma=m*n*(N^2-4)/(48*(N-1));
else
   mod(N,2)==1;
    mu=n*(N^2-1)/(4*N);
    sigma=m*n*(N+1)*(N^2+3)/(48*N^2);
end
t2k=(t2k-mu)^2/sigma;
end