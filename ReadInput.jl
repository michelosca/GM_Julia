module ReadInput

using SharedData
using GenerateODEs
using SharedData: me, mp, kb, e
using SharedData: Reaction, Species
using SharedData: AddSpeciesToList, AddReactionToList
using SharedData: SetSpeciesID
using SharedData: s_electron_id, s_Ar_id, s_ArIon_id

using PowerInput: PowerInputFunction

# Constants
const c_io_block_end = 0::Int64
const c_io_error = 1::Int64
const c_io_block_start = 4::Int64
const c_io_entry = 2::Int64
const c_io_emptyline = 3::Int64

function SetupInputData()

    errcode = c_io_error

    errcode = ReadInputData("input.deck")
    if (errcode != 0)
        print("***ERROR*** Failed to read the input deck. Abort code\n")
        return errcode
    end

    return errcode
end

global block_id = 0
global b_system = 1
global b_species = 2
global b_reactions = 3

###############################################################################
# START BLOCKS
function ReadInputData(filename)

    errcode = 1

    # Opens the file given in filename and reads each line
    # First read
    print("First reading of the input deck...\n")
    global read_step = 1 
    global input_id = 0
    errcode = ReadFile(filename)
    if (errcode == c_io_error)
        return errcode
    end
    print("End of input deck reading\n\n")

    # Second read
    print("Second reading of the input deck...\n")
    global read_step = 2 
    global input_id = 0
    errcode = ReadFile(filename)
    if (errcode == c_io_error)
        return errcode
    end
    print("End of input deck reading\n\n")
    return errcode
end


function ReadFile(filename::String)
    
    errcode = 0

    open(filename,"r") do f
        line = 1
        while ! eof(f)
            s = readline(f)
            # ReadLine identifies each line on filename
            errcode = ReadLine(s)
            if (errcode == c_io_error)
                print("***ERROR*** Stop reading at file line ", line,"\n")
                return errcode  
            end
            line += 1
        end
    end
    return errcode
end


function ReadLine(string)

    errcode = 1

    # Trimm out commets
    i_comment = findfirst("#", string)
    if !(i_comment === nothing)
        i_comment = i_comment[1]
        string = string[begin:i_comment-1]
    end
    string = strip(string)

    # Signal wether this is a begin/end:block 
    i_block = findfirst(":", string)
    # Signla whether it is a "name = var" line
    i_eq = findfirst("=", string)
    # Signla whether it is a reaction line
    i_react = findfirst(";", string)

    # Check line
    if !(i_block === nothing)
        i_block = i_block[1]
        global block_name = string[i_block+1:end]
        if (occursin("begin", string))
            errcode = StartBlock(block_name)
            if (errcode == c_io_error)
                print("***WARNING***\n  - Something went wrong starting the ",
                    block_name, " block\n")
            end
        elseif (occursin("end", string))
            errcode = EndBlock(block_name)
            if (errcode == c_io_error)
                print("***WARNING***\n  - Something went wrong ending the ",
                    block_name, " block\n")
            end
        end
    elseif !(i_eq === nothing)
        i_eq = i_eq[1]
        name = strip(string[begin:i_eq-1])
        var = strip(string[i_eq+1:end])
        errcode = ReadInputDeckEntry(name, var, block_id)
        if (errcode == c_io_error)
            print("***WARNING***\n  - Entry in ", block_name,
                "-block has not been located\n")
            print("  - Input entry: ", name ," = ",var ,"\n")
        end
    elseif !(i_react === nothing)
        errcode = ReadInputDeckEntry("",string, block_id) 
        if (errcode == c_io_error)
            print("***WARNING***\n  - Entry in ", block_name,
                "-block has not been located\n")
            print("  - Input entry: ", string ,"\n")
        end
    else
        errcode = c_io_emptyline
    end
    return errcode 
end

###############################################################################
# BLOCK START/END
function StartBlock(name)

    errcode = c_io_error
    if (occursin("system",name))
        global block_id = b_system
        errcode = StartSystemBlock()
    elseif (occursin("species",name))
        global block_id = b_species
        errcode = StartSpeciesBlock()
    elseif (occursin("reactions",name))
        global block_id = b_reactions
        errcode = StartReactionsBlock()
    end
    return errcode
end

function EndBlock(name)

    errcode = c_io_error
    global block_id = 0
    if (occursin("system",name))
        errcode = EndSystemBlock()
    elseif (occursin("species",name))
        errcode = EndSpeciesBlock()
    elseif (occursin("reactions",name))
        errcode = EndReactionsBlock()
    end
    return errcode
end
#
###############################################################################
###############################################################################
# START BLOCKS
function StartSystemBlock()
    errcode = 0
    return errcode
end


function StartSpeciesBlock()

    errcode = 0
    global input_id += 1
    if (read_step == 2)
        global input_n_id= 0 
        global input_elastic_id= 0
        global input_mass = 0.0
        global input_charge = 0.0
        global input_neq_flag = false
        global input_Teq_flag = false
        global input_wl_flag = false
        global input_P_flag = false
    end
    return errcode
end

function StartReactionsBlock()
    errcode = c_io_error
    if (read_step == 1)
        global input_r_id = 0
        global input_n_id = 0
        global input_inv_s = Int64[]
        global input_bal_s = Int64[]
        global input_rea_s = Int64[]
        global input_rate_coeff = Expr 
        global input_E = 0.0
        if (SharedData.s_electron_id != 0)
            errcode = 0
        else
            print("***ERROR*** Electron species must be defined before reactions\n")
        end
    else
        errcode = 0
    end
    return errcode
end
#
###############################################################################
###############################################################################
# END BLOCKS
function EndSystemBlock()
    errcode = 0
    return errcode
end

function EndSpeciesBlock()

    errcode = c_io_error 

    if (read_step == 1)
        errcode = 0
    elseif (read_step == 2)
        if (input_charge != 0)
            if (input_Teq_flag)
                global input_P_flag = true
            end
            global input_wl_flag = true
        end

        errcode = SharedData.AddSpeciesToList(input_id, input_mass,
            input_charge, input_neq_flag, input_Teq_flag, input_wl_flag,
            input_P_flag, input_n_id, input_elastic_id)
    end
    return errcode 
end

function EndReactionsBlock()
    errcode = 0 
    return errcode
end
#
###############################################################################
###############################################################################
# READ INPUT ENTRIES 
function ReadInputDeckEntry(name, var, block_id)

    errcode = 1 
    
    if (block_id == b_system)
        errcode = ReadSystemEntry(name, var)
    elseif (block_id == b_species)
        errcode = ReadSpeciesEntry(name, var)
    elseif (block_id == b_reactions)
        errcode = ReadReactionsEntry(name, var)
    end

    return errcode
end

###############
### SYTEM BLOCK
function ReadSystemEntry(name, var)

    errcode = 1

    if (read_step == 1)
        units_fact = 1.0
        units_index = findlast("_", var)
        if (units_index === nothing)
            units = ""
        else
            units_index = units_index[1]
            units_str = lowcase(var[units_index+1:end])
            if (units_str=="kw" || units_str=="khz")
                units_fact = 1.e3
                var = var[1:units_index-1]
            elseif (units_str=="mhz" || units_str=="mw")
                units_fact = 1.e6
                var = var[1:units_index-1]
            elseif (units_str=="ghz")
                units_fact = 1.e9
                var = var[1:units_index-1]
            end
        end
        
        errcode = SharedData.SetSystemParameters(name, var, units_fact)
    else
        errcode = 0
    end
    return errcode 
end

##################
### SPECIES BLOCK 
function ReadSpeciesEntry(name, var)

    errcode = c_io_error 

    if (name=="name")
        if (read_step == 1)
            errcode = SetSpeciesID(var, input_id)
        else 
            if (var=="e" || var=="electrons" || var=="electron")
                global input_n_id = s_electron_id
                errcode = 0
            elseif (occursin("Ar",var))
                global input_n_id = s_Ar_id
                errcode = 0
            else
                print("***ERROR*** Neutral species id has not been found\n")
                errcode = c_io_error 
            end
        end
    end

    if (name=="charge")
        if (read_step == 2)
            global input_charge = parse(Int64,var) * e
        end
        errcode = 0
    end

    if (name=="mass")
        if (read_step == 2)
            expr = Meta.parse(var)
            global input_mass = eval(expr)
        end
        errcode = 0
    end

    if (name=="solve_dens")
        if (read_step == 2)
            global input_neq_flag = parse(Bool, var) 
        end
        errcode = 0
    end

    if (name=="solve_temp")
        if (read_step == 2)
            global input_Teq_flag = parse(Bool, var)
        end
        errcode = 0
    end

    if (name=="dominant_reaction")
        if (read_step == 2)
            global input_elastic_id = parse(Int64, var)
        end
        errcode = 0
    end
    return errcode 
end

#################
### REACTION BLOCK
function ReadReactionsEntry(name, var)
    # Splits the input line contained in var into four parts
    # The fourth part is actually not necessary, but it is 
    # recommended to be included

    errcode = 0 

    if (read_step == 1)
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

        errcode = ParseReaction(reaction_process_str)
        if (errcode == c_io_error) return errcode end
        
        errcode = ParseEThreshold(e_threshold_str)
        if (errcode == c_io_error) return errcode end
        
        errcode = ParseRateCoefficient(rate_coeff_str)
        if (errcode == c_io_error) return errcode end

        if !(description_str === nothing)
            errcode = ParseDescription(description_str)
            if (errcode == c_io_error) return errcode end
        end

        errcode = SharedData.AddReactionToList(input_r_id, input_inv_s,
            input_rea_s, input_bal_s, input_rate_coeff, input_E,
            input_n_id)
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
        s_id = SharedData.s_Ar_id
        global input_n_id = s_id
    elseif (s == "Ar+")
        s_id = SharedData.s_ArIon_id
    elseif (s == "Ar*")
        s_id = SharedData.s_ArExc_id
    elseif (s == "e")
        s_id = SharedData.s_electron_id
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
        units_fact = 1
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
        global input_e_threshold = parse(Float64, str) * units_fact
    catch
        print("***ERROR*** While parsing E threshold\n")
        errcode = c_io_error
    end
        
    return errcode
end


function ParseRateCoefficient(str)

    # Must provide: rate_coefficient
    errcode = 0

    try
        expr = Meta.parse(str)
        errcode = SetRateCoefficientExpr(expr)
    catch
        errcode = c_io_error
        print("***ERROR*** While parsing the rate coefficient\n")
    end

    return errcode
end

function SetRateCoefficientExpr(expr)

    errcode = 0

    symbol_list = Symbol[]
    AssessExpression!(symbol_list, expr)
    AssessSymbols!(symbol_list)

    # Identify the Symbols involved
    symbol_found = false
    find_Te = findall(x-> x==:Te, symbol_list)
    if length(find_Te) > 0
        errcode = GenerateODEs.SetSymbolFlag(1)
        if (errcode==1) return errcode end
    end

    find_Te_eV = findall(x-> x==:Te_eV, symbol_list)
    if length(find_Te_eV) > 0
        errcode = GenerateODEs.SetSymbolFlag(2)
        if (errcode==1) return errcode end
    end

    find_m_Ar = findall(x-> x==:m_Ar, symbol_list)
    if length(find_m_Ar) > 0
        errcode = GenerateODEs.SetSymbolFlag(3)
        if (errcode==1) return errcode end
    end

    global input_rate_coeff = expr
    return errcode
end


function AssessExpression!(symbol_list, expr)
    # Identifies all the symbols in the given expression

    for e in expr.args
        if typeof(e) == Symbol
            # Bedore adding e to the list, check that it is not
            # already in symbol_list
            already_included = false
            for s in symbol_list
                if s == e
                    already_included = true
                    break
                end
            end
            if !already_included
                push!(symbol_list, e)
            end
        elseif typeof(e) == Expr
            AssessExpression!(symbol_list, e)
        end
    end
end


function AssessSymbols!(symbol_list)

    new_symbol_list = Symbol[]
    for s in symbol_list
        if s == :Te
            push!(new_symbol_list, s)
        elseif s == :Te_eV
            push!(new_symbol_list, s)
        elseif s == :m_Ar
            push!(new_symbol_list, s)
        end
    end
    symbol_list = new_symbol_list
end


function ParseDescription(str)

    # Must provide: id
    errcode = c_io_error

    str = lowercase(str)

    if (read_step == 1)
        if (str == "elastic")
            global input_r_id = SharedData.r_elastic_id 
            errcode = 0
        elseif (str == "excitation")
            global input_r_id = SharedData.r_exictat_id
            errcode = 0
        elseif (str == "ionization" || str == "ionisation")
            global input_r_id = SharedData.r_ionizat_id
            errcode = 0
        elseif (str == "recombination")
            global input_r_id = SharedData.r_recombi_id
            errcode = 0
        elseif (str == "charge exchange" || str == "charge-exchange" ||
            str == "cx" || str == "charge_exchange")
            global input_r_id = SharedData.r_cx_id
            errcode = 0
        end
    else
        errcode = 0
    end
    return errcode
end

end