# bond type for visualization
@data BondTypeVis begin
    Single
    Double
    Triple
    Aromatic
end

function _bond_marker(btype::BondTypeVis)
    @match btype begin
        Single => Single_Bond
        Double => Double_Bond
        Triple => Triple_Bond
        Aromatic => Aromatic_Bond
    end
    #return Single_Bond
end

function _bond_type(s::AbstractSystem, atom_id1::Integer, atom_id2::Integer)
    bo = bond_order(s, atom_id1, atom_id2)
    return if bo == 2//1
        Double
    elseif bo == 3//1
        Triple
    elseif bo == 3//2
        Aromatic
    else
        Single
    end
end

# data structure for bond position, direction, color, and type

Base.@kwdef mutable struct BondData
    origin::Vector{Point3f0} = Point3f[]
    direction::Vector{Point3f0} = Point3f[]
    colors::Vector{HSLA{Float32}} = HSLA{Float32}[]
end

function add_bdata!(b::BondData, origin::AbstractVector{F}, direction::AbstractVector{F}, color::HSLA{Float32}) where {F<:AbstractFloat}
    push!(b.origin, origin)
    push!(b.direction, direction)
    push!(b.colors, color)
    return nothing
end

function bondscatter!(axis::LScene, bonds::BondData, bond_radius::Number, marker::AbstractMesh)
    # scaling with absolute size
    scales = Vec3f[]
    for v in bonds.direction
        l = norm(v)
        push!(scales, Makie.Vec3f(bond_radius, bond_radius, l))
    end

    # rotation for bond
    rots = normalize(bonds.direction)

    meshscatter!(
        axis, bonds.origin;
        color = bonds.colors, marker = marker,
        rotation = rots, markersize = scales
    )
end

function bondscatter!(axis::LScene, bonds::Observable{BondData}, bond_radius::Number, marker::AbstractMesh)
    # scaling with absolute size
    scales = lift(bonds) do stub
        map(bonds[].direction) do v
            l = norm(v)
            Makie.Vec3f(bond_radius, bond_radius, l)
        end
    end

    # bond direction
    rots = lift(bonds) do stub
        normalize(bonds[].direction)
    end

    points = lift(bonds) do stub
        bonds[].origin
    end

    colors = lift(bonds) do stub
        bonds[].colors
    end

    meshscatter!(
        axis, points;
        color = colors, marker = marker,
        rotation = rots, markersize = scales
    )
end

function bond_nonpbc(s::AbstractSystem, color_func::Function)
    if wrapped(s)
        error("bond_nonpbc does not support wrapped system. Use bond_pbc instead.")
    end

    bonds = Dict{BondTypeVis, BondData}()
    bonds[Single] = BondData()
    bonds[Double] = BondData()
    bonds[Triple] = BondData()
    bonds[Aromatic] = BondData()
    for edge in edges(topology(s))
        atom_id1, atom_id2 = src(edge), dst(edge)
        p1, p2 = position(s, atom_id1), position(s, atom_id2)
        center = p1 + 0.5*(p2 - p1)
        bond_type = _bond_type(s, atom_id1, atom_id2)
        add_bdata!(
            bonds[bond_type],
            p1, center-p1, color_func(s, atom_id1)
        )
        add_bdata!(
            bonds[bond_type],
            center, p2-center, color_func(s, atom_id2)
        )
    end

    return bonds
end

function bond_pbc(s::AbstractSystem, color_func::Function)
    if !wrapped(s)
        error("bond_pbc does not support non-wrapped system. Use bond_nonpbc instead.")
    end

    # bond data structure
    bonds = Dict{BondTypeVis, BondData}()
    bonds[Single] = BondData()
    bonds[Double] = BondData()
    bonds[Triple] = BondData()
    bonds[Aromatic] = BondData()

    origin = box(s).origin # box origin
    axis = box(s).axis
    a, b, c = axis[:, 1], axis[:, 2], axis[:, 3] # box vector
    for edge in edges(topology(s))
        atom_id1, atom_id2 = src(edge), dst(edge)
        p1, p2 = position(s, atom_id1), position(s, atom_id2)
        diff = travel(s, atom_id2) - travel(s, atom_id1)
        bond_vector = (diff[1]*a + diff[2]*b + diff[3]*c + p2) - p1
        @assert all(x -> x ∈ (-1, 0, 1), diff)

        # solve p1 + t(p2 - p1) == αa + βb + γc + origin
        # ボンドが周期境界をまたぐときα, β, γのいずれかが1となる
        # ↓progress: t, cross_dim: 周期強化に交差する次元(1,2,3)，交差しなければ0
        cross_points = NamedTuple{(:progress, :cross_dim), Tuple{Float64, Int64}}[]
        resize!(cross_points, count(!=(0), diff)+3)
        #bondの始点，中点，終点
        cross_points[1] = (progress=0.0, cross_dim=0)
        cross_points[end-1] = (progress=0.5, cross_dim=0)
        cross_points[end] = (progress=1.0, cross_dim=0)
        n = 2
        if diff[1] != 0
            # x軸方向の周期境界とbondが交差
            α = diff[1]
            (_t, β, γ) = (-p1 + α*a + origin) \ [bond_vector -b -c] # wrong
            cross_points[n] = (progress=_t, cross_dim=1)
            n += 1
        end
        if diff[2] != 0
            β = diff[2]
            (_t, α, γ) = (-p1 + β*b + origin) \ [bond_vector -a -c]
            cross_points = (progress=_t, cross_dim=2)
            n += 1
        end
        if diff[3] != 0
            γ = diff[3]
            (_t, α, β) = (-p1 + γ*c + origin) \ [bond_vector -a -b]
            cross_points = (progress=_t, cross_dim=3)
            n += 1
        end
        sort!(cross_points; by=x->x.progress)
        @assert all(x -> -1 <= x.progress <= 1, cross_points)

        println(cross_points)

        # 結合次数別にbonds::Vector{BondData}を構築し，bond dataを追加
        color1, color2 = color_func(s, atom_id1), color_func(s, atom_id2)
        bond_type = _bond_type(s, atom_id1, atom_id2)
        pbc_shift = zeros(3)
        for i in 1:length(cross_points)-1
            t, t_next = cross_points[i].progress, cross_points[i+1].progress
            cl = t < 0.5 ? color1 : color2
            p  = p1 + t_next*bond_vector + pbc_shift
            println(p)
            add_bdata!(
                bonds[bond_type],
                p, (t-t_next)*bond_vector, cl
            )
            dim = cross_points[i+1].cross_dim
            if dim != 0
                pbc_shift -= diff[dim] * axis[:, dim]
            end
        end
    end

    return bonds
end
