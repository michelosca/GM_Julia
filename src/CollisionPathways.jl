module CollisionPathways

using Plots
using Plots.PlotMeasures
#using LaTeXStrings
using SharedData: Species, OutputBlock, Reaction, SpeciesID, System
using SharedData: K_to_eV, kb
using PlotModule: FindReactingSpecies, FindProductSpecies, Get2DMap
using DataFrames, CSV
using OutputModule: copy_species
using Printf


mutable struct CollisionPathway
    name::String
    pathway::Vector{Species}
    asc_prob::Float64
    des_prob::Float64
    
    CollisionPathway() = new()
end 

function InitializeCollisionPathway()
    new_pathway = CollisionPathway()
    new_pathway.name = ""
    new_pathway.pathway = Species[]
    new_pathway.asc_prob = 1.0
    new_pathway.des_prob = 1.0
    return new_pathway
end

function CopyCollisionPathway(cp::CollisionPathway)
    new_cp = InitializeCollisionPathway()
    new_cp.name = cp.name
    for node in cp.pathway
        push!(new_cp.pathway, copy_species(node))
    end
    new_cp.asc_prob = copy(cp.asc_prob)
    new_cp.des_prob = copy(cp.des_prob)
    return new_cp
end

function MainCollisionPathways(species_list::Vector{Species},
    reaction_list::Vector{Reaction}, sID::SpeciesID, source::String,
    data_param::Vector{Tuple{Symbol, Float64}},
    start_end_species_ids::Vector{Int64}, node_list::Vector{Int64})

    # Get node links between species
    GenerateNodeLinks!(species_list, reaction_list, sID, source, data_param)

    # Get collision pathways
    #start_end_species_ids = [sID.O2, sID.O_3s]
    #node_list = [sID.O2, sID.O, sID.O_1d, sID.O_3s, sID.O_3p, sID.O_5p, sID.O_5s]
    pathway_list = GetCollisionPathways(start_end_species_ids, node_list,
        species_list)

    # Get total probabilities for each pathway
    for pathway in pathway_list
        GetPathwayProbability!(pathway)
    end

    # Sort pathway list using the descending probability (top to bottom)
    SortCollisionPathways!(pathway_list)

    for pathway in pathway_list
        if pathway.des_prob >= 0.0001
            @printf("%100s; P_des[*100] = %6.4f\n", pathway.name, pathway.des_prob*100.0)
        end
    end
    return pathway_list
end

function SortCollisionPathways!(pathway_list::Vector{CollisionPathway})

    # Sort from max.-to-min. des.prob.
    len_list = length(pathway_list)
    sort_not_finished = true
    while sort_not_finished
        sort_not_finished = false
        for i in range(2,len_list,step=1)
            prev_pw = CopyCollisionPathway(pathway_list[i-1])
            pw = CopyCollisionPathway(pathway_list[i])
            #print("Prev prob ", prev_pw.des_prob)
            #print(" prob ", pw.des_prob, "\n")
            if pw.des_prob > prev_pw.des_prob
                #print("Switch\n")
                pathway_list[i-1] = pw
                pathway_list[i] = prev_pw
                sort_not_finished = true
            end
        end
    end

    #for pw in pathway_list
    #    @printf("%50s with des. prob. %06.4f\n", pw.name, pw.des_prob)
    #end
end

function GetPathwayProbability!(current_pathway::CollisionPathway)

    #print(current_pathway.name,"\n")
    
    ############################################################################
    # Calculate ascending probability
    current_pathway.asc_prob = 1.0 
    ix_path = 2
    for node in current_pathway.pathway[1:end-1]
        next_node_tuple = (-1,"None",0.0)
        for next_node in node.out_nodes
            next_node_id = next_node[1]
            next_species_id = current_pathway.pathway[ix_path].id
            if next_node_id == next_species_id 
                next_node_tuple = next_node
                break
            end
        end

        # Update next node index in pathway list
        ix_path += 1
        
        # Add node jump probability to total probability
        jump_probability = next_node_tuple[3]
        current_pathway.asc_prob *= jump_probability 
        #@printf("  - From %9s to %9s with ascending prob. of %04.2f\n", node.name, i_tuple[2], i_tuple[3])
    end
    #@printf("  - TOTAL ascending probability of %06.4f\n", current_pathway.asc_prob)

    ############################################################################
    # Calculate desceinding probability
    current_pathway.des_prob = 1.0 
    ix_path = length(current_pathway.pathway) - 1
    for node in Iterators.reverse(current_pathway.pathway[2:end])
        prev_node_tuple = (-1,"None",0.0)
        for prev_node in node.in_nodes
            prev_node_id = prev_node[1]
            prev_species_id = current_pathway.pathway[ix_path].id
            if prev_node_id == prev_species_id 
                prev_node_tuple = prev_node
                break
            end
        end

        # Update next node index in pathway list
        ix_path -= 1
        
        # Add node jump probability to total probability
        jump_probability = prev_node_tuple[3]
        current_pathway.des_prob *= jump_probability 
        #@printf("  - From %9s to %9s with ascending prob. of %04.2f\n", node.name, i_tuple[2], i_tuple[3])
    end
    #@printf("  - TOTAL descending probability of %06.4f\n", current_pathway.des_prob)

end

function GenerateNodeLinks!(species_list::Vector{Species},
    reaction_list::Vector{Reaction}, sID::SpeciesID, source::String,
    data_param::Vector{Tuple{Symbol, Float64}})

    prob_threshold = 0.0001

    # Load rate coefficient data frame
    K_df = CSV.read(source * "K_vs_P_Ar_vs_P%_O2_vs_power.csv", DataFrame)
    K_df_subset = DF_filter(K_df, data_param)

    for s1 in species_list
        # Reset in/out node lists
        s1.in_nodes = Tuple{Int64, String, Float64}[]
        s1.out_nodes = Tuple{Int64, String, Float64}[]

        if s1.id == sID.electron
            continue
        end
        # Find reactions that have the species in the LHS 
        s1_reactions_LHS = FindReactingSpecies(reaction_list, s1)
        # Find reactions that have the species in the RHS 
        s1_reactions_RHS = FindProductSpecies(reaction_list, s1)

        # Get total reaction rates
        K1_total_LHS = 0.0 #zeros(nrows, ncols)
        for r in s1_reactions_LHS 
            K_name = "r" * string(r.id) * ": "*r.name
            K1_total_LHS += K_df_subset[1, K_name]
        end
        K1_total_RHS = 0.0 #zeros(nrows, ncols)
        for r in s1_reactions_RHS 
            K_name = "r" * string(r.id) * ": "*r.name
            K1_total_RHS += K_df_subset[1, K_name]
        end

        # Loop over all species and find the in/out-node species
        for s2 in species_list
            if s2.id == sID.electron
                continue
            #elseif s2.id == s1.id
            #    continue
            end

            # Ascending probability: rate of s2 that is produced by s1
            s2_reactions_RHS = FindProductSpecies(s1_reactions_LHS, s2)
            K2_total_RHS = 0.0 #zeros(nrows, ncols)
            for r in s2_reactions_RHS
                K_name = "r" * string(r.id) * ": "*r.name
                K2_total_RHS += K_df_subset[1, K_name]
            end
            ascend_rate = K2_total_RHS / K1_total_LHS 
            if ascend_rate >= prob_threshold
                push!(s1.out_nodes, (s2.id, s2.name, ascend_rate))
            end

            # Descending probability: rate of s1 that is produced by s2
            s2_reactions_LHS = FindReactingSpecies(s1_reactions_RHS, s2)
            K2_total_LHS = 0.0 #zeros(nrows, ncols)
            for r in s2_reactions_LHS
                K_name = "r" * string(r.id) * ": "*r.name
                K2_total_LHS += K_df_subset[1, K_name]
            end
            descend_rate = K2_total_LHS / K1_total_RHS 
            if descend_rate >= prob_threshold
                push!(s1.in_nodes, (s2.id, s2.name, descend_rate))
            end
        end 

        # Print results 
        #print(s1.name,"\n")
        #for node in s1.in_nodes
        #    @printf(" - Incoming %9s with probability %05.3f\n",node[2], node[3])
        #end
        #for node in s1.out_nodes
        #    @printf(" - Outcoming %9s with probability %05.3f\n",node[2], node[3])
        #end
    end
end

function GetCollisionPathways(start_end_species::Vector{Int64},
    node_list::Vector{Int64}, main_species_list::Vector{Species})

    # Set start/end species nodes
    start_species_id = start_end_species[1]
    end_species_id = start_end_species[2]
    start_species = main_species_list[start_species_id]
    end_species = main_species_list[end_species_id]

    # Get node-set 
    species_list = Species[]
    for s in main_species_list
        if s.id == start_species.id
            continue
        elseif s.id == end_species.id
            continue
        end
        for s_id in node_list
            if s_id == s.id
                push!(species_list, s)
                break
            end
        end
    end
    nodes = length(species_list)

    # This is the output argument: a list of pathway lists
    pathway_list = CollisionPathway[]

    # This function outputs a list with the number of pathways sorted by number of nodes
    # For instance, a set of 2 nodes has the following combinations:
    # - 1 pathway with zero nodes
    # - 2 pathways with 1 nodes
    # - 2 pathways with 2 nodes
    # This is gathered in a list that is pathway_options = [1,2,2]
    pathway_options = GetMaxNumberPathways(nodes)

    # Loop over the pathway n-node options
    node_jumps = 0
    for n_node_paths in pathway_options
        #print(" - Pathways with ", node_jumps, " node jumps\n")

        # path_track keeps track of the different path combinations that may be possible
        path_track = ones(Int64, node_jumps) 

        # Loop over the number of pathways there are with n_node_step node-jumps
        for i_path in range(1,n_node_paths, step=1)

            # Parameters used in case a node reaches a dead-end
            valid_path = true
            last_jump = node_jumps

            # Array where the species pathway is saved
            current_pathway = InitializeCollisionPathway()

            # First node in the pathway is the start species
            push!(current_pathway.pathway, start_species)
        
            # The pathway must have "node_jumps" steps
            for jump in range(1, node_jumps, step=1)

                # The next-node options come from the previous node
                current_node = current_pathway.pathway[end]
                
                # Loops over the nex-node options until a valid new node is found
                next_node_not_found = true
                while next_node_not_found
                    #print("    Path track ", path_track)
                    #print("; Leading species ", current_node.name)

                    # If path_track at a given node > than number of available node-options
                    # then this pathway has reached a dead-end
                    if path_track[jump] > length(current_node.out_nodes)
                        #print("\nDead end path\n")
                        valid_path = false
                        last_jump = jump
                        break
                    end
                    next_tuple = current_node.out_nodes[path_track[jump]] 
                    next_speciesnode = main_species_list[next_tuple[1]]
                    #print("; new node species ", next_speciesnode.name,"\n")

                    # Discard species if:
                    #  - end-species
                    end_species_node_flag = next_speciesnode.id == end_species_id
                    #  - start-species
                    start_species_node_flag = next_speciesnode.id == start_species_id
                    # - Already in the current pathway 
                    node_in_currentpath = FindNodeInNodeset(current_pathway.pathway, next_speciesnode)
                    # - Species is no part of the node set under consideration 
                    node_in_nodelist = FindNodeInNodeset(species_list, next_speciesnode)
                    if end_species_node_flag || start_species_node_flag ||
                        node_in_currentpath || !node_in_nodelist
                        #print("Species discarded\n")
                        # If node is discarded then move to next available node
                        path_track[jump] += 1 
                        continue
                    end
                    
                    # Available node is a good option: add to the current pathway and exit while loop
                    push!(current_pathway.pathway, next_speciesnode)
                    next_node_not_found = false
                end # while loop

                # If pathway option reached a dead-end then exit the current_pathway loop
                if !valid_path
                    break
                end

            end # jump loop


            # Update step
            # If it was a dead-end, set last path_track index to 1
            if !valid_path
                if last_jump > 1
                    path_track[last_jump] = 1
                end
                last_jump -= 1
            end

            # Update path_track
            for i in last_jump:-1:1 
                path_track[i] += 1
                pathway_ix = i + 1
                len_current_pathway = length(current_pathway.pathway[pathway_ix-1].out_nodes)
                path_track_num = path_track[i]
                if path_track_num > len_current_pathway 
                    if i != 1
                        path_track[i] = 1
                    end
                else
                    break
                end
            end

            # If dead-end found at first path_track index then this pathway completelly discarded
            if !valid_path 
                if path_track[1] > length(start_species.out_nodes)
                    break
                else
                    continue
                end
            end

            # Pathway is valid if end-species is an option as a next node
            #print("Last chaing species ", current_pathway[end].name, "\n")
            valid_end_species = FindNextNode(current_pathway.pathway[end], end_species)
            if valid_end_species
                #print("End species added\n")
                push!(current_pathway.pathway, end_species)
            else
                #print("End species NOT added\n")
                valid_path = false 
            end
            
            # Print pathway 
            if valid_path
                pathway_str = GeneratePathwayString(current_pathway.pathway)
                current_pathway.name = pathway_str
                push!(pathway_list, current_pathway)
                #print(current_pathway.name, "\n")
            end

            if node_jumps > 0
                # Stop trying more possible pathways if first index path_track is larger
                # than the node-options in the start species
                if path_track[1] > length(start_species.out_nodes)
                    break
                end
            end
        end
        node_jumps += 1

    end

    return pathway_list
end

function GeneratePathwayString(pathway::Vector{Species})
    pathway_str = ""
    for node in pathway
        pathway_str *= @sprintf("%s -> ", node.name)
    end
    pathway_str = pathway_str[1:end-4]
    return pathway_str
end

function GetMaxNumberPathways(nodes::Int64)

    pathway_list = zeros(Int64, nodes + 1)

    ix = 1
    pathway_list[ix] = 1

    sum = 0
    for i in range(1,nodes,step=1)
        prod = 1
        for j in range(1,i-1,step=1)
            prod *= nodes - j
        end

        ix += 1
        pathway_list[ix] = prod * nodes

        sum += prod
    end
    return pathway_list 
end

function FindNextNode(species::Species, next_species::Species)

    new_species_id = next_species.id
    for new_node in species.out_nodes
        new_node_id = new_node[1]
        if new_species_id == new_node_id 
            return true
        end
    end
    return false
end

function FindNodeInNodeset(species_list::Vector{Species}, species::Species)

    species_id = species.id
    for s in species_list
        if s.id == species_id
            return true
        end
    end
    return false
end

function DF_filter(df::DataFrame, data_param::Vector{Tuple{Symbol, Float64}})
    
    df_subset = df
    n_cols = length(data_param)
    for c in range(1,n_cols, step=1)
        data_tuple = data_param[c]
        col = data_tuple[1] 
        val = data_tuple[2] 
        df_data = FindMatchingDFValue(df, col, val)
        df_subset = df_subset[ df_subset[!,col] .== df_data, :]
    end
    return df_subset 
end

function FindMatchingDFValue(df::DataFrame, col_name::Symbol, value::Union{Float64,Int64})
    data = df[!,col_name]
    len_data = length(data)

    matching_value = 0.0
    match_found = false
    
    i = 1
    min_val = data[i]
    max_val = 0.5 * (data[i+1] + data[i])
    if min_val <= value && value < max_val
        matching_value = min_val
        return matching_value
    end

    for i in range(2,len_data-1,step=1)
        current_data = data[i]
        min_val = 0.5 * (data[i-1] + current_data)
        max_val = 0.5 * (data[i+1] + current_data)
        if min_val <= value && value < max_val
            matching_value = current_data 
            match_found = true 
            break
        end
    end

    if !match_found
        i = len_data 
        min_val = 0.5 * (data[i-1] + data[i])
        max_val = data[i]
        if min_val <= value && value <= max_val
            matching_value = max_val
            match_found = true
        end
    end
    if !match_found
        print("*** WARNING *** Match value was not found\n")
    end
    return matching_value
end

end