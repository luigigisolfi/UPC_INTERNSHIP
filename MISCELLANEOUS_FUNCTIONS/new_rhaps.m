function NR = new_rhaps(func, param)
    % Call fzero with the function handle and the initial guess separately
    NR = fzero(func, param); 
end


