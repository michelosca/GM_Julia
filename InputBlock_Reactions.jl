module InputBlock_Reactions

using SharedData: c_io_error, e, K_to_eV
using SharedData: Reaction
using InputBlock_Species: s_electron_id
using InputBlock_Species: s_Ar_id, s_ArIon_id, s_ArExc_id
using InputBlock_Species: s_O_id, s_O2_id, s_OIon_id, s_OnIon_id, s_O2Ion_id
using InputBlock_Species: s_O3p_id, s_O1d_id, s_O2a1Ag_id
try
    using ReactionSet: K_funct_list
    print("ReactionSet loaded\n")
catch
    print("ReactionSet not yet loaded\n")
end

###############################################################################
################################  VARIABLES  ##################################
###############################################################################
global reaction_list = Reaction[]

# REACTION IDs 
# In case more reaction need to be defined, they need to be added here
const r_case_energy_sink = -1
const r_wall_loss = -2
const r_elastic_id = 1
const r_excitat_id = 2
const r_ionizat_id = 3
const r_recombi_id = 4
const r_cx_id = 5

###############################################################################
################################  FUNCTIONS  ##################################
###############################################################################
# FUNCTION TREE
# - StartFile_Reactions
# - StartReactionsBlock
# - ReadReactionsEntry
#   - ParseReaction
#     - GetSpeciesFromString
#       - SelectSpeciesID
#     - GetReactionSpeciesLists
#   - ParseEThreshold
#   - ParseRateCoefficient
#     - ReplaceSymbol
#   - ParseDescription
#   - AddReactionToList 
# - EndReactionsBlock

function StartFile_Reactions(read_step) 

    errcode = 0

    if (read_step == 1)
        global reaction_list = Reaction[]
    end

    return errcode
end

function StartReactionsBlock(read_step)
    errcode = c_io_error

    global r_count = 0
    if (read_step == 1)
        global input_r_id = 0
        global input_n_id = 0
        global input_inv_s = Int64[] 
        global input_bal_s = Int64[]
        global input_rea_s = Int64[]
        global input_E = 0.0
        global input_r_case = 0
        if (s_electron_id != 0)
            errcode = 0
        else
            print("***ERROR*** Electron species must be defined before reactions\n")
        end
    elseif (read_step == 0)
        global f_ReactionSet = open("ReactionSet.jl","w") 
        open("ReactionSet.Template","r") do f_temp
            while ! eof(f_temp)
                line_str = readline(f_temp, keep = true)
                if (line_str == "### START REACTION STRINGS ###\n")
                    print("Found start reaction string\n")
                    break
                else
                    write(f_ReactionSet,line_str)
                end
            end
        end
        errcode = 0
    else
        errcode = 0
    end
    return errcode
end


function ReadReactionsEntry(name, var, read_step)
    # Splits the input line contained in var into four parts
    # The fourth part is actually not necessary, but it is 
    # recommended to be included

    errcode = 0 
    global r_count += 1

    if (read_step == 1 || read_step == 0)
        # First part: the reaction process
        idx = findfirst(";", var)
        if !(idx === nothing)
            idx = idx[1]
            reaction_process_str = strip(var[1:idx-1])
            var = var[idx+1:end]
            str_track = true
        else
            errcode = c_io_error
            str_track = false 
        end

        # Second part: the threshold energy
        idx = findfirst(";", var)
        if (!(idx === nothing) && str_track)
            idx = idx[1]
            e_threshold_str = strip(var[1:idx-1])
            var = var[idx+1:end]
            str_track = true
        else
            errcode = c_io_error
            str_track = false 
        end

        # Third part: the rate coefficient expression
        idx = findfirst(";", var)
        if (!(idx === nothing) && str_track)
            idx = idx[1]
            rate_coeff_str = strip(var[1:idx-1])
            var = var[idx+1:end]
            str_track = true
        else
            errcode = c_io_error
            str_track = false
        end

        # Fourth part, if existing: reaction description string
        if (str_track)
            description_str = strip(var)
        else
            description_str = nothing 
        end
        if (errcode == c_io_error) return errcode end

        # Parse each term of the input line
        if (read_step ==1)
            errcode = ParseReaction(reaction_process_str)
            if (errcode == c_io_error) return errcode end
            
            errcode = ParseEThreshold(e_threshold_str)
            if (errcode == c_io_error) return errcode end

            if !(description_str === nothing)
                errcode = ParseDescription(description_str)
                if (errcode == c_io_error) return errcode end
            end
            
            global input_K_funct = K_funct_list[r_count]
        
            errcode = AddReactionToList(input_r_id, input_inv_s,
                input_rea_s, input_bal_s, input_K_funct, input_E,
                input_n_id, input_r_case)
        elseif (read_step == 0)
            errcode = WriteRateCoefficientsToModule(rate_coeff_str)
            if (errcode == c_io_error) return errcode end
            end
    end

    return errcode
end


function ParseReaction(str)

    # Must provide: neutral_species_id, involved_species,
    #               species_balance, reactant_species
    errcode = c_io_error

    idx = findfirst("->",str)
    if !(idx===nothing)
        # Get reactant and product strings
        reactant_str = strip(str[1:idx[1]-1])
        product_str = strip(str[idx[2]+1:end])
        
        # Get Vector{Int64} with species IDs
        input_rea_s = GetSpeciesFromString(reactant_str)
        if (input_rea_s == c_io_error) return errcode end

        input_pro_s = GetSpeciesFromString(product_str)
        if (input_pro_s == c_io_error) return errcode end

        # Get balance, reactant and involved species vectors
        errcode = GetReactionSpeciesLists(input_rea_s, input_pro_s)
        if (input_pro_s == c_io_error) return errcode end
    end
    return errcode
end


function GetSpeciesFromString(str)
    # This function links the species found in the given string (str)
    # to the species ids predefined in the code

    s_list = Int64[]
    next_species = true 
    while next_species
        idx = findfirst(" + ", str)
        if (idx===nothing)
            s_id = SelectSpeciesID(str)
            if (s_id == 0)
                print("***ERROR*** Reaction species ",str ," is not recognized\n")
                return c_io_error
            else 
                push!(s_list, s_id)
            end
            next_species = false
        else
            s = strip(str[1:idx[1]-1])
            str = strip(str[idx[2]+1:end])
            s_id = SelectSpeciesID(s)
            if (s_id == 0)
                print("***ERROR*** Reaction species ",s ," is not recognized\n")
                return c_io_error
            else 
                push!(s_list, s_id)
            end
        end
    end
    return s_list
end


function SelectSpeciesID(s)

    s_id = 0
    if (s == "Ar")
        s_id = s_Ar_id
        global input_n_id = s_id
    elseif (s == "Ar+")
        s_id = s_ArIon_id
    elseif (s == "Ar*")
        s_id = s_ArExc_id
    elseif (s == "e")
        s_id = s_electron_id
    elseif (s == "O")
        s_id = s_O_id
        global input_n_id = s_id
    elseif (s == "O2")
        s_id = s_O2_id
        global input_n_id = s_id
    elseif (s == "O+")
        s_id = s_OIon_id
    elseif (s == "O-")
        s_id = s_OnIon_id
    elseif (s == "O2+")
        s_id = s_O2Ion_id
    elseif (s == "O(3p)")
        s_id = s_O3p_id
    elseif (s == "O(1d)")
        s_id = s_O1d_id
    elseif (s == "O2(a1Ag)")
        s_id = s_O2a1Ag_id
    end
    return s_id
end


function GetReactionSpeciesLists(reac::Vector{Int64}, prod::Vector{Int64})

    errcode = 0
    try
        involved_species = Int64[]
        reactant_species = Int64[]

        # Loop over species list and
        #  - if not in involved_species -> push
        #  - if already in, move on
        # Loop for reactant species
        for r in reac
            already_in = false
            i = 0
            for is in involved_species
                i += 1
                if r == is
                    already_in = true
                    break
                end
            end
            if !already_in
                push!(involved_species, r)
                push!(reactant_species, r)
            end
        end
        
        # Can add reactant species to global parameter
        global input_rea_s = reactant_species 

        # Loop for product species
        for p in prod 
            already_in = false
            i = 0
            for is in involved_species
                i += 1
                if p == is
                    already_in = true
                    break
                end
            end
            if !already_in
                push!(involved_species, p)
            end
        end

        # Can add involved species to global parameter
        global input_inv_s = involved_species

        # Get species balance between reactants and products
        n_species = length(involved_species)
        balance_list = zeros(Int64, n_species) 
        for i in 1:n_species
            s = involved_species[i]
            for r in reac
                if s==r
                    balance_list[i] -= 1
                end
            end

            for p in prod
                if s==p
                    balance_list[i] += 1
                end
            end
        end

        # Add balance list to global parameter
        global input_bal_s = balance_list 
    catch
        print("***ERROR*** While setting up reaction species lists\n")
        errcode = c_io_error
    end
    return errcode
end


function ParseEThreshold(str)

    # Must provide: E_threshold
    errcode = 0

    try
        units_fact = 1.0
        if (occursin("J", str))
            idx = findfirst("J", str)
            idx = idx[1]-1
            str = string(str[1:idx])
        elseif (occursin("eV", str))
            units_fact = e
            idx = findfirst("eV", str)
            idx = idx[1]-1
            str = string(str[1:idx])
        end
        global input_E = parse(Float64, str) * units_fact
    catch
        print("***ERROR*** While parsing E threshold\n")
        errcode = c_io_error
    end
        
    return errcode
end


function WriteRateCoefficientsToModule(str)

    # Must provide: rate_coefficient
    errcode = 0 

    try
        expr = Meta.parse(str)
        if !(typeof(expr)==Float64)
            ReplaceSymbolS!(expr)
        end

        # Now that the expression is ready to be evaluated, write it down in a new file
        if input_r_id == r_wall_loss
            write(f_ReactionSet, string("push!(K_funct_list, (temp, species) -> ",expr,")\n"))
        else
            write(f_ReactionSet, string("push!(K_funct_list, (temp) -> ",expr,")\n"))
        end

    catch
        errcode = c_io_error
        print("***ERROR*** While writing rate coefficient to module\n")
    end

    return errcode
end


function ReplaceSymbolS!(expr::Expr)
    # In case new symbols are used in the input.deck,
    # they would need to be added here
    ReplaceSymbol!(expr, :Te_eV,     Expr(:call,:*,:(temp[s_electron_id]), K_to_eV))
    ReplaceSymbol!(expr, :Te,        :(temp[s_electron_id]))
    ReplaceSymbol!(expr, :TArIon_eV, Expr(:call,:*,:(temp[s_ArIon_id]), K_to_eV))
    ReplaceSymbol!(expr, :TArIon, :(temp[s_ArIon_id]))
    ReplaceSymbol!(expr, :TO2,   :(temp[s_O2_id]))
    ReplaceSymbol!(expr, :TO,    :(temp[s_O_id]))
    ReplaceSymbol!(expr, :m_Ar,  Expr(:call,:*,:amu, 40 ))
    ReplaceSymbol!(expr, :m_O,   Expr(:call,:*,:amu, 16 ))
    ReplaceSymbol!(expr, :m_O2,  Expr(:call,:*,:amu, 32 ))
    ReplaceSymbol!(expr, :m_O2,  Expr(:call,:*,:amu, 32 ))
    ReplaceSymbol!(expr, :R, :(system.radius) )
    ReplaceSymbol!(expr, :L, :(system.l) )
    ReplaceSymbol!(expr, :A, :(system.A) )
    ReplaceSymbol!(expr, :V, :(system.V) )
    ReplaceSymbol!(expr, :h_R, :(species.h_R) )
    ReplaceSymbol!(expr, :h_L, :(species.h_L) )
    ReplaceSymbol!(expr, :uB, :(species.u_Bohm) )
    ReplaceSymbol!(expr, :vth, :(species.u_thermal) )
    ReplaceSymbol!(expr, :Lambda, :(species.Lambda) )
    ReplaceSymbol!(expr, :gamma, :(species.gamma) )
    ReplaceSymbol!(expr, :D, :(species.D) )
end


function ReplaceSymbol!(expression::Expr, old, new)
    n = length(expression.args)
    for i in 1:n
        arg = expression.args[i]
        if (typeof(arg) == Expr)
            ReplaceSymbol!(arg, old, new)
        elseif (arg == old)
            expression.args[i] = new
        end
    end
end


function ParseDescription(str)

    # Must provide? id
    errcode = 0 

    str = lowercase(str)

    if (str == "elastic")
        global input_r_id = r_elastic_id 
        errcode = 0
    elseif (str == "excitation")
        global input_r_id = r_excitat_id
        errcode = 0
    elseif (str == "ionization" || str == "ionisation")
        global input_r_id = r_ionizat_id
        errcode = 0
    elseif (str == "recombination")
        global input_r_id = r_recombi_id
        errcode = 0
    elseif (str == "charge exchange" || str == "charge-exchange" ||
        str == "cx" || str == "charge_exchange")
        global input_r_id = r_cx_id
        errcode = 0
    elseif (str == "excitation energy_sink" || str == "excitation energy sink")
        global input_r_id = r_excitat_id
        global input_r_case = r_case_energy_sink
    elseif (str == "wall_rate_coefficient")
        global input_r_id = r_wall_loss
        errcode = 0
    end

    return errcode
end


function AddReactionToList(id::Int64, invol_s::Vector{Int64},
    react_s::Vector{Int64}, balan_s::Vector{Int64}, K_funct, E::Float64,
    n_id::Int64, r_case::Int64)

    errcode = 0
    try
        # Add react to reaction_list
        react = Reaction(id, n_id, invol_s, balan_s, react_s, K_funct, E, r_case) 
        push!(reaction_list, react)
    catch
        print("***ERROR*** While attaching reaction\n")
        errcode = c_io_error 
    end
    return errcode
end


function EndReactionsBlock(read_step)
    errcode = 0 

    if (read_step == 0)
        write_flag = false 
        open("ReactionSet.Template","r") do f_temp
            while ! eof(f_temp)
                line_str = readline(f_temp, keep = true)

                if write_flag
                    write(f_ReactionSet,line_str)
                end

                if (line_str == "### END REACTION STRINGS ###\n")
                    print("Found end reaction string\n")
                    write_flag = true
                end
                    
            end
        end
        close(f_ReactionSet)
    end

    return errcode
end

end