module HMDVisualize

using Colors
# Choose colorschemes with care! Refer to Peter Kovesi's PerceptualColourMaps package, or to Fabio Crameri's Scientific Colour Maps for more information.
using ColorSchemes: colorschemes, get
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
using Symbolics

export visualize, color_scheme, viridis, default_color, atom_color

###
###### preset colors
###

const atom_color = Dict(
    elements[:H ].number => colorant"hsla(  0,   0%,  90%, 1.0)",
    elements[:C ].number => colorant"hsla(249,  14%,  50%, 1.0)",
    elements[:N ].number => colorant"hsla(240,  70%,  78%, 1.0)",
    elements[:O ].number => colorant"hsla(351,  95%,  70%, 1.0)",
    elements[:F ].number => colorant"hsla( 95,  40%,  50%, 1.0)",
    elements[:Si].number => colorant"hsla( 20,  39%,  74%, 1.0)",
    elements[:S ].number => colorant"hsla( 41,  95%,  59%, 1.0)",
    elements[:Cl].number => colorant"hsla(158,  79%,  42%, 1.0)",
    elements[:Br].number => colorant"hsla(  0,  60%,  36%, 1.0)"
)

function default_color(s::AbstractSystem, atom_id::Integer)
    elem = element(s, atom_id)
    return if elem >= 1
        atom_color[elem]
    else elem < 0
        colorant"hsla(170,  0%, 37%, 1.0)"
    end
end

include("assets/bonds.jl")
include("bond.jl")

###
###### system visualization functions
###

function visualize(s::AbstractSystem{D, F, SysType}; color_func::Function=default_color, atom_radius::Number=0.30, bond_radius::Number=0.15) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return visualize(
        Trajectory(s);
        color_func = color_func,
        atom_radius = atom_radius,
        bond_radius = bond_radius
    )
end

#TODO: 回転中心の指定，平行移動速度の自動調整，回転速度のスライドバー指定
function visualize(traj::AbstractTrajectory{D, F, SysType}; color_func::Function=default_color, atom_radius::Number=0.30, bond_radius::Number=0.15) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    if dimension(traj[1]) != 3
        error("expected dimension 3, found $D")
    end
    wrap_coord = wrapped(traj)

    # enumerate all possible colors
    #colors = Vector{HSLA{Float32}}(undef, 0)
    #for reader in traj
    #    snapshot = reader.reader
    #    for atom_id in 1:natom(snapshot)
    #        color = color_func(snapshot, atom_id)
    #        if color ∉ colors
    #            push!(colors, color)
    #        end
    #    end
    #end

    fig = Figure(; backgroundcolor = :black)
    sl_x = Slider(fig[2, 1], range = 1:length(traj), startvalue = 1)
    axis = LScene(fig[1,1]; show_axis = false)
    cam3d!(axis; projectiontype = :orthographic, mouse_translationspeed=0.001f0)

    reader = similar_system(traj)

    # box plot
    box_mesh = lift(sl_x.value) do index
        update_reader!(reader, traj, index)
        get_boxmesh(reader)
    end
    mesh!(axis, box_mesh; color = :white)

    # atom plot
    atoms = lift(box_mesh) do stub
        Point3f.(all_positions(reader))
    end
    atom_colors = lift(box_mesh) do stub
        [color_func(reader, atom_id) for atom_id in 1:natom(reader)]
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
    bonds = lift(box_mesh) do stub
        if wrapped(traj)
            bond_pbc(reader, color_func)
        else
            bond_nonpbc(reader, color_func)
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
    return HSLA{Float32}(get(colorschemes[scheme], value))
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
