module ReadInput

include("SharedData.jl")
using .SharedData

export GetInputData
export InitializeData

# Constants
const c_io_block_end = 0
const c_io_block_start = 1
const c_io_entry = 2
const c_io_emptyline = 3

function InitializeData()

    # Define species ids
    #c_electron_id = 1
    ar_ion_id = 3

    # Define reaction ids
    r_ela_id = 1
    r_ion_id = 2
    r_rec_id = 3

    # Define electron-argon elastic scattering
    f_r1(temp::Vector{Float64}) = 2.336e-14 * temp[c_electron_id]^1.609
    Er1(temp::Vector{Float64}) = 3 * me / mg * kb * (temp[c_electron_id] - temp[c_neutral_id])
    r1 = Reaction(r_ela_id, [c_electron_id,c_neutral_id],[c_electron_id, c_neutral_id],[0,0], f_r1, Er1)

    # Define e-impact ionization: e + Ar -> e + e + Ar+
    f_r3(temp::Vector{Float64}) = 2.34e-14 * temp[c_electron_id]^0.59 * exp(-17.44/temp[c_electron_id])
    Er3(temp::Vector{Float64}) =  15.76 * e
    r3 = Reaction(r_ion_id, [c_electron_id,c_neutral_id,ar_ion_id],[c_electron_id, c_neutral_id], [1,-1,1], f_r3, Er3)
    
    # Define Ar recombination: e + Ar+ -> Ar
    f_r4(temp::Vector{Float64}) =  5e-39 * temp[c_electron_id]^4.5
    Er4(temp::Vector{Float64}) =  0
    r4 = Reaction(r_rec_id, [c_electron_id,ar_ion_id, c_neutral_id], [c_electron_id, ar_ion_id], [-1,-1,1], f_r4, Er4)

    # Load species
    f_zero(dens::Vector{Float64} ,temp::Vector{Float64}) = 0
    fwl_n_e(dens::Vector{Float64} ,temp::Vector{Float64}) = 0
    fwl_T_e(dens::Vector{Float64} ,temp::Vector{Float64}) = 0
    fpi_n_e(dens::Vector{Float64} ,temp::Vector{Float64}) = 0
    fwl_n_i(dens::Vector{Float64} ,temp::Vector{Float64}) = 0
    fwl_T_i(dens::Vector{Float64} ,temp::Vector{Float64}) = 0
    fpi_n_i(dens::Vector{Float64} ,temp::Vector{Float64}) = 0
    #                   species_id   , n_eq,  T_eq,  wl_n  , fwl_n  , wl_T  , fwl_T  , pi_n ,  fpi_n
    electrons = Species(c_electron_id, true,  true , false , fwl_n_e, false , fwl_T_e, false,  fpi_n_e)
    argon =     Species(c_neutral_id , true,  false, false , f_zero , false , f_zero , false,  f_zero)
    argon_io =  Species(ar_ion_id    , true,  false, true  , fwl_n_i, true  , fwl_T_i, false,  f_zero)

    species_list = [ electrons, argon, argon_io ]
    reaction_list = [r1,r3,r4]

    return species_list, reaction_list

end

function GetInputData(filename)

    # Opens the file given in filename and reads each line
    try
        open(filename,"r") do f
            line = 0
            while ! eof(f)
                s = readline(f)
                # ReadLine identifies each line on filename
                read_flag = ReadLine(s)
                line += 1
            end
        end
        return 1
    catch
        print("Input deck does not exist!\n")
        return 0
    end

end

function ReadLine(string)

    # Trimm out commets
    i_comment = findfirst("#", string)
    if !(i_comment === nothing)
        i_comment = i_comment[1]
        string = string[begin:i_comment-1]
    end

    # Check wether this is a begin/end:block 
    i_block = findfirst(":", string)
    if !(i_block === nothing)
        i_block = i_block[1]
        block_name = string[i_block+1:end]
        if (occursin("begin", string))
            print(block_name, " ",c_io_block_start,"\n")
            return c_io_block_start
        elseif (occursin("end", string))
            print(block_name, " ", c_io_block_end,"\n")
            return c_io_block_end
        end
    else
        # Check whether it is a "name = var" line
        i_eq = findfirst("=", string)
        if !(i_eq === nothing)
            i_eq = i_eq[1]
            name = string[begin:i_eq-1]
            var = string[i_eq+1:end]
            print("Name: ", name, "\n")
            print(" Var: ", var , "\n\n")
            return c_io_entry
        end
    end
    return c_io_emptyline

end

#function ParseFunction(string)
#    expr = Meta.parse(string)
#    funct = @eval (param_list) -> $expr
#end
end