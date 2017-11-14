%% perform the outer product of 2 matrices
function C = outer_prod(A,B)
dimA = size(A);
dimB = size(B);
C = zeros(dimA(1)*dimB(1),dimA(2)*dimB(2));
for i = 1:dimA(1)
    for j = 1:dimA(2)
        C((i-1)*dimB(1)+1:i*dimB(1),(j-1)*dimB(2)+1:j*dimB(2)) = A(i,j)*B;
    end
end

