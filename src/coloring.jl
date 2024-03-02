function count_ltypes(s::AbstractSystem{D, F, S}) where {D, F<:AbstractFloat, S<:AbstractSystemType}
    count = counter(String)
    for label in all_labels(s, "polymeric")
        inc!(count, type(label))
    end
    return count
end


struct MoleculeViridis <: StaticColoring
    colors::Vector{LCHuvA{Float32}}
end
atom_colors(mv::MoleculeViridis) = mv.colors

function MoleculeViridis(s::AbstractSystem{D, F, S}) where {D, F<:AbstractFloat, S<:AbstractSystemType}
    mols = all_labels(s, "polymeric", Polymer)
    contiguous = all(i -> id(mols[i]) + 1 == id(mols[i+1]), 1:length(mols)-1)

    colors = Vector{LCHuvA{Float32}}(undef, natom(s))
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



struct VisualizeMolecule <: StaticColoring
    colors::Vector{LCHuvA{Float32}}
end
atom_colors(vm::VisualizeMolecule) = vm.colors

function VisualizeMolecule(s::System, molecule::HLabel)
    mol_atoms = sort!(label2atom(s, "polymeric", molecule))
    colors = map(1:natom(s)) do atom
        if insorted(atom, mol_atoms)
            return default_color(s, atom)
        else
            return palette[:transparent]
        end
    end
    return VisualizeMolecule(colors)
end


struct VisualizeRU <: StaticColoring
    colors::Vector{LCHuvA{Float32}}
end
atom_colors(vm::VisualizeRU) = vm.colors

function VisualizeRU(s::System, ru::HLabel)
    ru_atoms = sort!(label2atom(s, "polymeric", ru))
    polymer_atoms = let
        super_molecule = super(s, "polymeric", ru)
        sort!(label2atom(s, "polymeric", super_molecule))
    end

    colors = map(1:natom(s)) do atom
        if insorted(atom, ru_atoms)
            return default_color(s, atom)
        elseif insorted(atom, polymer_atoms)
            dc = default_color(s, atom)
            return LCHuvA{Float32}(dc.l, dc.c, dc.h, 0.3)
        else
            return palette[:transparent]
        end
    end
    return VisualizeRU(colors)
end


struct EmphAtoms <: StaticColoring
    colors::Vector{LCHuvA{Float32}}
end
atom_colors(ea::EmphAtoms) = ea.colors

function EmphAtoms(s::System, atoms::Union{Base.Generator, <:AbstractRange, Vector{<:Integer}})
    _atoms = sort(atoms)
    colors = map(1:natom(s)) do atom
        if insorted(atom, _atoms)
            return palette[:emph]
        else
            return default_color(s, atom)
        end
    end
    return EmphAtoms(colors)
end

function EmphAtoms(s::System, atoms::Union{Base.Generator, <:AbstractRange, Vector{<:Integer}}, label::HLabel)
    _atoms = sort(atoms)
    label_atoms = sort!(label2atom(s, "polymeric", label))
    colors = map(1:natom(s)) do atom
        if insorted(atom, _atoms)
            return palette[:emph]
        elseif insorted(atom, label_atoms)
            return default_color(s, atom)
        else
            return palette[:transparent]
        end
    end
    return EmphAtoms(colors)
end
