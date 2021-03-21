%Expands a contact matrix, so that age-stratified matrix is extended so
%that its aged groups are further sub-stratified uniformly. E.g. An age
%group of 0-1 years may be sub-stratified into 12 monthly groups
%0-1month, 1-2month, etc. It expands from the leftmost column up to
%divide column and from the topmost row to the divide row.
%A - Contact Matrix. Must be square. Age groups must be evenly spaced.
%In the from aij is the contacts of age group i with age group j.
%divide - Dividing column and row. Must be positive integer.
%expand - The amount that each age-group want to be sub-divided into. Must
%be positive integer. E.g. if 12 each age-group will become 12 age-groups.
function x = uniformsubstratify(A,divide,expand)
q1 = A(1:divide,1:divide);
q2 = A(1:divide,divide+1:end);
q3 = A(divide+1:end,1:divide);
q4 = A(divide+1:end,divide+1:end);
%Divide q1 and q3 by the expansion number, to around for uniformity of
%liklihood of contact amoungst age-groups.
q1 = q1./expand;
q3 = q3./expand;
%Transform each entry in q1 into a expand by expand ones matrix and combine
%into matrix in the correct order.
[q1length,q1width] = size(q1);
newq1 = [];
for i=1:q1length
    newrow=[];
    for j=1:q1width
        value = q1(i,j);
        newmatrix = value*ones(expand);
        newrow = [newrow,newmatrix];
    end
    newq1 = [newq1;newrow];
    newrow = [];
end
%Elongate each entry in q2 into a vector expand long and combine back into
%a matrix in the correct order.
[q2length,q2width] = size(q2);
newq2 = [];
for i=1:q2length
    newrow=[];
    for j=1:q2width
        value = q2(i,j);
        newmatrix = value*ones(expand,1);
        newrow = [newrow,newmatrix];
    end
    newq2 = [newq2;newrow];
    newrow = [];
end
%Widen each entry in q3 into row vector expand wide and combine back into a
%matrix in the correct order.
[q3length,q3width] = size(q3);
newq3 = [];
for i=1:q3length
    newrow=[];
    for j=1:q3width
        value = q3(i,j);
        newmatrix = value*ones(1,expand);
        newrow = [newrow,newmatrix];
    end
    newq3 = [newq3;newrow];
    newrow = [];
end
%Recombine the matrix.
new1 = [newq1,newq2];
new2 = [newq3,q4];
x = [new1;new2];