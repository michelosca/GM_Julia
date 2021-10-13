module ReadInput

export GetInputData

# Constants
const c_io_block_end = 0
const c_io_block_start = 1
const c_io_entry = 2
const c_io_emptyline = 3

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