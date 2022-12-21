function [i, j] = linear_sum_assignment(cost_matrix)

True = 1;
False = 0;

if length(size(cost_matrix)) ~= 2
    error(sprintf('Expected a matrix (2-d array), got a %i-d array'), ...
          length(size(cost_matrix))); 
end

%  The algorithm expects more columns than rows in the cost matrix.
if size(cost_matrix, 2) < size(cost_matrix, 1)
    cost_matrix = cost_matrix';
    transposed = True;
else
    transposed = False;
end

state = Hungary_(cost_matrix); 

step = @step1; 

while ~isequal(step, @none)
    [step, state] = step(state); 
end

if transposed
    marked = state.marked';
else
    marked = state.marked;
end
    
[i, j] = find(marked == 1); 


function s = Hungary_(cost_matrix)
[n, m] = size(cost_matrix); 

s = struct('C', cost_matrix, ...
           'row_uncovered', true(n, 1), ...
           'col_uncovered', true(1, m), ...
           'Z0_r', 1, ... 
           'Z0_c', 1, ... 
           'path', ones(n+m, 2, 'int64'), ...
           'marked', zeros(n, m, 'int8')); 

function none()
return 

function state = clear_covers(state)
True = 1;
False = 0;

state.row_uncovered(:) = True; 
state.col_uncovered(:) = True;



function [fnhandle, state] = step1(state)
True = 1;
False = 0; 

state.C = state.C - ...
          repmat(min(state.C,[],2), [1, size(state.C, 2)]); 

[i_, j_] = find(state.C == 0);

for k = 1:length(i_)
    i = i_(k); 
    j = j_(k); 
    if state.col_uncovered(j) & state.row_uncovered(i)
        state.marked(i, j) = 1;
        state.col_uncovered(j) = False;
        state.row_uncovered(i) = False;
    end
end

state = clear_covers(state);
fnhandle = @step3;

function [fnhandle, state] = step3(state)
True = 1;
False = 0;
marked = (state.marked == 1);
state.col_uncovered( any(marked) ) = False; 

if sum(marked(:)) < size(state.C, 1)
    fnhandle = @step4; 
else
    fnhandle = @none; 
end
        
function [fnhandle, state] = step4(state)
True = 1;
False = 0;

[n, m] = size(state.C); 

covered_C = logical(state.C == 0);
covered_C(~state.row_uncovered, :) = 0;
covered_C(:, ~state.col_uncovered, :) = 0;

while True
    idx_zero = find(covered_C, 1); 
    
    if length(idx_zero) == 0
        fnhandle = @step6; 
        return 
    else
        row = mod(idx_zero-1, n)+1; col = floor((idx_zero-1)/n)+1; 
        state.marked(row, col) = 2;
        %# Find the first starred element in the row
        star_col = find(state.marked(row, :) == 1, 1); 
        if isempty(star_col)
            %# Could not find one
            state.Z0_r = row;
            state.Z0_c = col;
            fnhandle = @step5; 
            return
        else
            col = star_col;
            state.row_uncovered(row) = False;
            state.col_uncovered(col) = True;
            covered_C(~state.row_uncovered, col) = 0 ; 
            covered_C(row, :) = 0; 
        end
    end
end


function [fnhandle, state] = step5(state)
True = 1;
False = 0;
el0 = 1;
el1 = 2;
[n, m] = size(state.C); 

count = 1; 
state.path(count, el0) = state.Z0_r;
state.path(count, el1) = state.Z0_c;

while True
    %# Find the first starred element in the col defined by
    %# the path.
    row = find(state.marked(:, state.path(count, el1) ) == 1, 1); 
    if isempty(row)
        %# Could not find one
        break; 
    else
        count = count + 1; 
        state.path(count, el0) = row; 
        state.path(count, el1) = state.path(count - 1, el1);
    end
        
    col = find(state.marked( state.path(count, el0), : ) == 2, 1); 
    count = count + 1; 
    state.path(count, el0) = state.path(count - 1, el0);
    state.path(count, el1) = col;
end
    
for i = 1:(count)
    if state.marked( state.path(i, el0), state.path(i, el1) ) == 1
        state.marked( state.path(i, el0), state.path(i, el1)) = 0; 
    else
        state.marked( state.path(i, el0), state.path(i, el1) ) = 1;
    end
end
    
state = clear_covers(state);
%# Erase all prime markings
state.marked(state.marked == 2) = 0;
fnhandle = @step3; return 


function [fnhandle, state] = step6(state)
if any(state.row_uncovered) & any(state.col_uncovered)
    minval = min(min(state.C(state.row_uncovered, ...
                             state.col_uncovered))); 
    state.C(~state.row_uncovered, :) =  state.C(~state.row_uncovered, :) ...
        + minval;
    state.C(:, state.col_uncovered) = state.C(:, state.col_uncovered) ...
        - minval; 
end

fnhandle = @step4;