n = 6;
attempts = 0;
while(true)
    A = rand(n,n) > 0.7;
    A = triu(A, 1) + ones(size(A, 1)).*triu(A, 1)'
    G = graph(A);
    if ~all(conncomp(G) == 1)
        display('not connected');
        bins = conncomp(G)
        for i=1:length(bins)
            if bins(i) ~= 1
                A(i, i-1) = 1;
                A(i-1, i) = 1;
                G = graph(A);
                bins = conncomp(G)
            end
        end
    end
    plot(G)
    if all(1)%cycleCountBacktrack('adjMatrix', A) == 0)
        display(['Number of failed attempts', num2str(attempts)])
        break
    else
        attempts = attempts + 1;
    end
end