module HMDVisualize

using Colors
# Choose colorschemes with care! Refer to Peter Kovesi's PerceptualColourMaps package, or to Fabio Crameri's Scientific Colour Maps for more information.
using ColorSchemes: colorschemes, get
using DataStructures
using GeometryBasics
using Graphs
using HMD
using HMDPolymer
using LinearAlgebra
using LoopVectorization
using MLStyle
using GLMakie #cite
using PeriodicTable
using Printf
using StaticArrays
using Symbolics

export visualize, color_scheme, viridis, atom_color
export AbstractAtomColoring, StaticColoring, DefaultColoring, MoleculeViridis, count_ltypes
export VisualizeMolecule, VisualizeRU, EmphAtoms

###
###### preset colors
###

const atom_color = Dict(
    elements[:H ].number => LCHuvA{Float32}(91.1, 0.00 , 140 , 1.0),
    elements[:C ].number => LCHuvA{Float32}(47.9, 28.4 , 271 , 1.0),
    elements[:N ].number => LCHuvA{Float32}(68.4, 65.6 , 266 , 1.0),
    elements[:O ].number => LCHuvA{Float32}(63.6, 104  , 5.57, 1.0),
    elements[:F ].number => LCHuvA{Float32}(66.9, 67.6 , 116 , 1.0),
    elements[:Si].number => LCHuvA{Float32}(75.9, 28.5 , 37.2, 1.0),
    elements[:P].number  => LCHuvA{Float32}(75.9, 28.5 , 37.2, 1.0), # temporary!
    elements[:S ].number => LCHuvA{Float32}(79.6, 92.0 , 55.7, 1.0),
    elements[:Cl].number => LCHuvA{Float32}(69.0, 68.7 , 148 , 1.0),
    elements[:Br].number => LCHuvA{Float32}(33.2, 82.6 , 12.2, 1.0)
)

const palette = Dict(
    :emph => LCHuvA{Float32}(69.6, 114  , 25.57, 1.0),
    :transparent => LCHuvA{Float32}(0, 0, 0, 0.0),
    :nonatom => LCHuvA{Float32}(40, 0, 214, 1.0)
)

function default_color(s::AbstractSystem, atom_id::Integer)
    elem = element(s, atom_id)
    return if elem >= 1
        atom_color[elem]
    else elem < 0
        palette[:nonatom]
    end
end

abstract type AbstractAtomColoring end
abstract type StaticColoring <: AbstractAtomColoring end

function (coloring::StaticColoring)(s::AbstractSystem{D, F, S}) where {D, F<:AbstractFloat, S<:AbstractSystemType}
    return atom_colors(coloring)
end

struct DefaultColoring{T<:AbstractSystemType} <: StaticColoring
    atom_color::Vector{LCHuvA{Float32}}
end
atom_colors(dc::DefaultColoring) = dc.atom_color

function DefaultColoring(
    s::AbstractSystem{D, F, S}
) where {D, F<:AbstractFloat, S<:AbstractSystemType}
    return DefaultColoring{S}([default_color(s, i) for i in 1:natom(s)])
end

function DefaultColoring(
    s::AbstractSystem{D, F, BeadsSpring}
) where {D, F<:AbstractFloat}
    nelem = maximum(all_elements(s))
    colors = [viridis(element(s, atom_id) / nelem) for atom_id in 1:natom(s)]
    return DefaultColoring{BeadsSpring}(colors)
end

include("coloring.jl")
include("assets/bonds.jl")
include("bond.jl")

###
###### system visualization functions
###

function visualize(
    s::AbstractSystem{D, F, SysType},
    fig = Figure(; backgroundcolor = :black),
    time_obs = nothing;
    atom_radius::Real = 0.30,
    bond_radius::Real = 0.15,
    boxvisualize::Bool = true,
    coloring::AbstractAtomColoring = DefaultColoring(s)
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    traj = Trajectory(s)
    return visualize(
        traj,
        fig,
        time_obs;
        atom_radius = atom_radius,
        bond_radius = bond_radius,
        boxvisualize = boxvisualize,
        coloring = coloring,
    )
end

# traits for BS system
struct NonBeadsSpring <: AbstractSystemType end
_is_BS(traj::Trajectory{D, F, BeadsSpring, L}) where {D, F, L} = BeadsSpring()
_is_BS(traj::Trajectory{D, F, <:AbstractSystemType, L}) where {D, F, L} = NonBeadsSpring()

function visualize(
    traj::AbstractTrajectory,
    fig = Figure(; backgroundcolor = :black),
    time_obs = nothing;
    atom_radius::Real = 0.30,
    bond_radius::Real = 0.15,
    boxvisualize::Bool = true,
    coloring::AbstractAtomColoring = DefaultColoring(traj[1])
)
    return visualize(
        _is_BS(traj),
        traj,
        fig,
        time_obs,
        atom_radius,
        bond_radius,
        boxvisualize,
        coloring
    )
end

# needs holy trait
function visualize(
    ::BeadsSpring,
    traj::Trajectory{D, F, BeadsSpring, L},
    fig,
    time_obs,
    atom_radius::Real,
    bond_radius::Real,
    boxvisualize::Bool,
    coloring::AbstractAtomColoring
) where {D, F<:AbstractFloat, L}
    return _visualize(traj, fig, time_obs, atom_radius, bond_radius, boxvisualize, coloring)
end

function visualize(
    ::NonBeadsSpring,
    traj::Trajectory{D, F, SysType, L},
    fig::Makie.Figure,
    time_obs,
    atom_radius::Real,
    bond_radius::Real,
    boxvisualize::Bool,
    coloring::AbstractAtomColoring
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType, L}
    return _visualize(
        traj,
        fig,
        time_obs,
        atom_radius,
        bond_radius,
        boxvisualize,
        coloring
    )
end

#TODO: 回転中心の指定，平行移動速度の自動調整，回転速度のスライドバー指定
function _visualize(
    traj::AbstractTrajectory{D, F, SysType},
    fig::Makie.Figure,
    time_obs,
    atom_radius::Real,
    bond_radius::Real,
    boxvisualize::Bool,
    coloring::AbstractAtomColoring
) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    if dimension(traj[1]) != 3
        error("expected dimension 3, found $D")
    end
    wrap_coord = wrapped(traj)

    #fig = Figure(; backgroundcolor = :black)
    #sl_x = Slider(fig[2, 1], range = 1:length(traj), startvalue = 1)
    if isnothing(time_obs)
        time_obs = Slider(fig[2, 1], range = 1:length(traj), startvalue = 1).value
    end
    axis = LScene(fig[1,1]; show_axis = false)
    cam3d!(axis; projectiontype = :orthographic, mouse_translationspeed=0.001f0, cad=false)

    reader = similar_system(traj)

    # box plot
    #box_mesh = lift(sl_x.value) do index
    box_mesh = lift(time_obs) do index
        update_reader!(reader, traj, index)
        get_boxmesh(reader)
    end
    mesh!(
        axis, box_mesh;
        color = boxvisualize ? :white : :transparent
    )

    # atom plot
    atoms = lift(box_mesh) do stub
        Point3f.(all_positions(reader))
    end
    atom_colors = lift(box_mesh) do stub
        colors = coloring(reader)
    end
    inspect_labels = map(1:natom(reader)) do i
        p = position(reader, i)
        elem = element(reader, i)
        elem = if elem in 1:length(elements)
            string(elements[elem].symbol)
        else
            string(elem)
        end
        x = @sprintf("%6f", p[1])
        y = @sprintf("%6f", p[2])
        z = @sprintf("%6f", p[3])
        "id: $i\n\
        element: $(elem)\n\
        x: $(x)\n\
        y: $(y)\n\
        z: $(z)"
    end
    meshscatter!(axis, atoms;
        color = atom_colors,
        markersize = atom_radius,
        inspector_label = (self, i, p) -> inspect_labels[i]
    )

    # bonds plot
    bonds = lift(atom_colors) do colors
        if wrapped(traj)
            bond_pbc(reader, colors)
        else
            bond_nonpbc(reader, colors)
        end
    end
    for bondtype in keys(bonds[])
        bonds_typed = @lift $(bonds)[bondtype]
        bondscatter!(axis, bonds_typed, bond_radius, _bond_marker(bondtype))
    end

    return fig
end

function get_boxmesh(s)
    a, b, c = box(s).axis[:,1], box(s).axis[:,2], box(s).axis[:,3]
    s_origin = box(s).origin
    p(i, j, k) = Point{3, Float32}(s_origin .+ (i .* a) .+ (j .* b) .+ (k .* c))

    #line_bewteen!(axis, p(0,0,0), p(0,1,0))
    #line_bewteen!(axis, p(0,0,0), p(1,0,0))
    #line_bewteen!(axis, p(0,1,0), p(1,1,0))
    #line_bewteen!(axis, p(1,0,0), p(1,1,0))
    #line_bewteen!(axis, p(0,0,1), p(0,1,1))
    #line_bewteen!(axis, p(0,0,1), p(1,0,1))
    #line_bewteen!(axis, p(0,1,1), p(1,1,1))
    #line_bewteen!(axis, p(1,0,1), p(1,1,1))
    #line_bewteen!(axis, p(0,0,0), p(0,0,1))
    #line_bewteen!(axis, p(1,0,0), p(1,0,1))
    #line_bewteen!(axis, p(0,1,0), p(0,1,1))
    #line_bewteen!(axis, p(1,1,0), p(1,1,1))

    return merge(normal_mesh.([
        Cylinder(p(0,0,0), p(0,1,0), Float32(0.06)),
        Cylinder(p(0,0,0), p(1,0,0), Float32(0.06)),
        Cylinder(p(0,1,0), p(1,1,0), Float32(0.06)),
        Cylinder(p(1,0,0), p(1,1,0), Float32(0.06)),
        Cylinder(p(0,0,1), p(0,1,1), Float32(0.06)),
        Cylinder(p(0,0,1), p(1,0,1), Float32(0.06)),
        Cylinder(p(0,1,1), p(1,1,1), Float32(0.06)),
        Cylinder(p(1,0,1), p(1,1,1), Float32(0.06)),
        Cylinder(p(0,0,0), p(0,0,1), Float32(0.06)),
        Cylinder(p(1,0,0), p(1,0,1), Float32(0.06)),
        Cylinder(p(0,1,0), p(0,1,1), Float32(0.06)),
        Cylinder(p(1,1,0), p(1,1,1), Float32(0.06))
    ]))
end

function line_bewteen!(axis, p1, p2)
    lines!(axis, [[p1[i], p2[i]] for i in 1:3]...)
    return nothing
end

function color_scheme(value::Real; scheme=:viridis)
    return LCHuvA{Float32}(get(colorschemes[scheme], value))
end

viridis(value::Real) = color_scheme(value::Real; scheme=:viridis)

function update_reader!(reader, traj, index)
    import_dynamic!(reader, traj, index)
    if is_reaction(traj, index)
        import_static!(reader, traj, index)
    end

    return nothing
end


end # module
