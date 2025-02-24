function name_value = struct_to_nv_pairs(s)
    % Helper function to convert a struct into name-value pairs.
    names = fieldnames(s);
    values = struct2cell(s);
    name_value = reshape([names, values]', 1, []);
end