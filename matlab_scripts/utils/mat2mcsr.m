%convert a matrix to modified compresses sparse row format
function y = mat2mcsr(A)

n = size(A, 1);

if size(A, 2) ~= n
    error('The matrix need to be square')
end

idx = zeros(3*n, 1);
val = zeros(3*n, 1);

%firstly write the diagonal part of the matrix
for i = 1:n
    val(i) = A(i, i);
end

%unsused value
val(n + 1) = 0;

cnt = n + 2;
%secondly traverse to off-diagonal element and write to idx and val
for i = 1:n
    idx(i) = cnt;
    for j = 1:n
        if j ~= i && A(i, j) ~= 0
            val(cnt) = A(i, j);
            idx(cnt) = j;
            cnt = cnt + 1;
        end
    end
end

idx(n + 1) = cnt;

y(1,:) = idx(1:cnt - 1) - 1;
y(2,:) = val(1:cnt - 1);
