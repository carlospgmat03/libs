__precompile__() 

module customs

using LinearAlgebra

export custom_unique

# Custom function to check if two complex numbers are approximately equal
function is_approx_equal(x, y)
    return norm(x - y) < tolerance
end

@doc "Customized unique function that aims to improve Julia's default unique function, a tolerance is used"
function custom_unique(list; tolerance=1e-8)
    filtered_vals = []
    for v in list
        if all(!is_approx_equal(v, fv) for fv in filtered_vals)
            push!(filtered_vals, v)
        end
    end
    return filtered_vals
end

end