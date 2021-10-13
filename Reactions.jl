module Reactions

function get_species_ids_and_count(s_list)
    # Input: s_list, is a species id list 
    # Output: list of species and their appearance

    involved_species_list = zeros(0) # Empty array tracking the involved species
    species_list = Tuple{Int, Int}[]
    for s in s_list 
        species_count = length(findall(x->x==s, s_list))
        if (length( findall( x -> x == s, involved_species_list )) == 0)
            append!(involved_species_list, s)
            push!(species_list, (s, species_count) )
        end
    end

    return species_list
end

end