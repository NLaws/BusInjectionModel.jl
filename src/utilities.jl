

function admittance_builder(net::Network{MultiPhase})

    order = bfs_order(net)
    Y = zeros(ComplexF64, (3*length(order), 3*length(order)))

    for edge in CommonOPF.edges(net)

        z_prim = CommonOPF.zij(edge[1],edge[2],net)

        phases = [1,2,3]
        missing_phases = []
        for i in [1,2,3]
            if z_prim[i,:] == ComplexF64[0.0 + 0.0im for i = 1:3]                
                filter!(!=(i),phases)
                append!(missing_phases,i)
            end
        end
        
        if length(phases) == 2
            z_prim = z_prim[1:end .!= missing_phases[1], 1:end .!= missing_phases[1]]
        elseif  length(phases) == 1
            z_prim = z_prim[1:end .== phases[1], 1:end .== phases[1]]
        end

        y_prim = inv(z_prim)

        position1 = findfirst(x -> x == edge[1],order)
        position2 = findfirst(x -> x == edge[2],order)

        for (row,phs1) in enumerate(phases)
            for (col,phs2) in enumerate(phases)
                Y[3*(position1-1)+phs1, 3*(position2-1)+phs2] += y_prim[row,col]
                Y[3*(position2-1)+phs2, 3*(position1-1)+phs1] += y_prim[row,col]
            end
        end
    end

    non_zero_rows = any(Y .!= 0.0 + 0.0im, dims=2)
    non_zero_indices = findall(non_zero_rows) .|> x -> x[1]
    Y = Y[non_zero_indices, non_zero_indices]

    return Y

end

function bfs_order(net::Network{MultiPhase})
    b = Graphs.bfs_tree(net.graph, 1)
    order = Any[1]
    i = 1
    while i <= length(order)
        append!(order, neighbors(b, order[i]))
        i += 1
    end

    for (i,num) in enumerate(order)
        order[i] = net.graph.vertex_labels[num]
    end

    return order
end