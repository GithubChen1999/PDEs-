
# f_gamete=["WW","WR","SW","SR"]
# m_gamete=["ww","wr","sw","sr"]



function allel2genotype(gamete)
    genotypes = [sort([f1, f2]) for f1 in gamete, f2 in gamete]  # Create all possible female genotypes by pairing gametes
    genotypes_unique = unique([g[1]* g[2] for g in genotypes])  # Convert the pairs to genotype strings and remove duplicates
    return genotypes_unique
end
f_genotypes_unique=allel2genotype(f_gamete)
m_genotypes_unique=allel2genotype(m_gamete)


function offspring_genotypes_number(f_gamete,m_gamete)    #generate the offspring poosible genotypes 
    offspring_genotype_number=[]
    for f in f_gamete, m in m_gamete
        f_s,m_s=sort([lowercase(f),m])
        push!(offspring_genotype_number, f_s*m_s)
    end
    return unique!(offspring_genotype_number)
end
offspring_genotype_number=offspring_genotypes_number(f_gamete,m_gamete)

function parents_g2_offspring_g(f_genotypes_unique,m_genotypes_unique)
    locus_number=Int8(length(f_genotypes_unique[1])/2)
    parents_g=[]
    for f_g in f_genotypes_unique, m_g in m_genotypes_unique
        push!(parents_g,[lowercase(f_g),m_g])
    end
    offspring_g=[]
    for idx in eachindex(parents_g)
        m_pair = [parents_g[idx][1][1:locus_number], parents_g[idx][1][locus_number+1:end]]   #male gametes
        f_pair = [parents_g[idx][2][1:locus_number], parents_g[idx][2][locus_number+1:end]]   #female gametes
        combined_pairs = [
        m_pair[1] * f_pair[1],  # ABab
        m_pair[1] * f_pair[2],  # ABcd
        m_pair[2] * f_pair[1],  # CDac
        m_pair[2] * f_pair[2]   # CDcd
        ]
        push!(offspring_g, combined_pairs)
    end
    for i in eachindex(offspring_g)
        for j in eachindex(offspring_g[i])
            pair=sort([offspring_g[i][j][1:locus_number],offspring_g[i][j][locus_number+1:end]]) 
            offspring_g[i][j]=pair[1]*pair[2]
        end
    end
    return offspring_g ,parents_g   
end
offspring_g, parents_g=parents_g2_offspring_g(f_genotypes_unique,m_genotypes_unique)


function genotype_count(offspring_genotype_number, offspring_g)
    genotype_count_matrix=fill(0,(length(offspring_genotype_number),length(offspring_g)))
    # Fill the matrix with counts of each genotype in the 100-element vector
    for (i, subvec) in enumerate(offspring_g)
        for (j, genotype) in enumerate(offspring_genotype_number)
            genotype_count_matrix[j, i] = count(x -> x == genotype, subvec)
        end
    end
    return genotype_count_matrix
end
genotype_count_matrix=genotype_count(offspring_genotype_number,offspring_g)

# display(offspring_genotype_number)
# display(parents_g)
# display(genotype_count_matrix)

total_unique_genotypes=unique(vcat(offspring_genotype_number,reduce(vcat,parents_g)))     # all genotype we get parents and offspring
if (length(total_unique_genotypes)-length(offspring_genotype_number)) >0
    gap=fill(0,(length(total_unique_genotypes)-length(offspring_genotype_number),length(total_unique_genotypes)))
    genotype_count_matrix=vcat(genotype_count_matrix,gap)
end
display(total_unique_genotypes)  #
display(parents_g)

m_prob_genotype=copy(genotype_count_matrix./8)  #y axis
f_prob_genotype=copy(genotype_count_matrix./8)   # xaxis

prob_genotype=copy(f_prob_genotype)
display(prob_genotype)
#for the x-shredder, the 
if x_shredder
    m_prob_genotype[:,2:3] .+=f_prob_genotype[:,2:3]
    f_prob_genotype[:,2:3].=0
    display(m_prob_genotype)
    display(f_prob_genotype)
end

#for YLE:
if YLE
    f_prob_genotype[:,2:3].=0
    m_prob_genotype[2,2] += m_prob_genotype[1,2]
    m_prob_genotype[1,2]=0
    display(m_prob_genotype)
    display(f_prob_genotype)
end

#for the fs-ridl gene:
 if fs_gene
    indices_with_t = Int8.(findall(x -> occursin("t", x), offspring_genotype_number))
    f_prob_genotype[indices_with_t,:].=0
    display(m_prob_genotype)
    display(f_prob_genotype)

end
#for early:
 J_fitness_matrix=sum(m_prob_genotype.+f_prob_genotype,dims=1)

