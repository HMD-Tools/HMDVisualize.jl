struct MoleculeViridis <: AbstractAtomColoring
    colors::Vector{LCHuv{Float32}}
end

function MoleculeViridis(s::AbstractSystem{D, F, S}) where {D, F<:AbstractFloat, S<:AbstractSystemType}
    mols = all_labels(s, "polymeric", Polymer)
    contiguous = all(i -> id(mols[i]) + 1 == id(mols[i+1]), 1:length(mols)-1)

    colors = Vector{LCHuv{Float32}}(undef, natom(s))
    for (i, mol) in enumerate(mols)
        atom_ids = label2atom(s, "polymeric", mol)
        if contiguous
            colors[atom_ids] = [viridis(id(mol) / length(mols)) for _ in atom_ids]
        else
            colors[atom_ids] = [viridis(i / length(mols)) for _ in atom_ids]
        end
    end

    return MoleculeViridis(colors)
end

function (coloring::MoleculeViridis)(s::AbstractSystem{D, F, S}) where {D, F<:AbstractFloat, S<:AbstractSystemType}
    return coloring.colors
end

function count_ltypes(s::AbstractSystem{D, F, S}) where {D, F<:AbstractFloat, S<:AbstractSystemType}
    count = counter(String)
    for label in all_labels(s, "polymeric")
        inc!(count, type(label))
    end
    return count
end
