function xop = read_ttop(d,file1)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
x = load(file1);
x = struct2cell(x);
i=1;
xop=tt_tensor(x{i},1e-14,size(x{i},2),size(x{i},1),size(x{i},3));
for i=2:d
    h=tt_tensor(x{i},1e-14,size(x{i},2),size(x{i},1),size(x{i},3));
    xop = tkron(xop,h);
end
%xop = tt_matrix(xop);
end
