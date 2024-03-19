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
Base.@kwdef struct BondData
    origin::Observable{Vector{Point3f0}} = Observable(Point3f[])
    direction::Observable{Vector{Makie.Vec3f}} = Observable(Makie.Vec3f[])
    scale::Observable{Vector{Makie.Vec3f}} = Observable(Makie.Vec3f[]) # (radius, radisu, length)
    colors::Observable{Vector{LCHuvA{Float32}}} = Observable(LCHuvA{Float32}[])
end

function BondData(n::Integer, bond_radius::Real)
    @assert n ≥ 0
    return BondData(
        Observable([Point3f(0.0f0, 0.0f0, 0.0f0) for _ in 1:n]),
        Observable([Makie.Vec3f(0.0f0, 0.0f0, 0.0f0) for _ in 1:n]),
        Observable([Makie.Vec3f(bond_radius, bond_radius, 0.0f0) for _ in 1:n]),
        Observable([palette[:emph] for _ in 1:n])
    )
end

function bonddata_by_type(s::AbstractSystem, bond_radius::Real)
    num_bond = Dict{BondTypeVis, Int64}(btype => 0 for btype in (Single, Double, Triple, Aromatic))
    for b in bonds(s)
        btype = _bond_type(s, src(b), dst(b))
        num_bond[btype] += 1
    end

    # bond multipliction for coloring and PBC
    n = wrapped(s) ? 5 : 2
    return Dict{BondTypeVis, BondData}(
        Single   => BondData(n * num_bond[Single]  , bond_radius),
        Double   => BondData(n * num_bond[Double]  , bond_radius),
        Triple   => BondData(n * num_bond[Triple]  , bond_radius),
        Aromatic => BondData(n * num_bond[Aromatic], bond_radius)
    )
end

function bondscatter!(axis::LScene, bonds::BondData, marker::AbstractMesh)
    meshscatter!(
        axis, bonds.origin;
        color = bonds.colors, marker = marker,
        rotation = bonds.direction, markersize = bonds.scale,
        inspectable = false
    )
end

function bond_nonpbc!(
    bdata::Dict{BondTypeVis, BondData},
    s::AbstractSystem{D, F, S},
    colors::AbstractVector{T},
    bondtypelist::Vector{BondTypeVis}
) where {D, F<:AbstractFloat, S<:AbstractSystemType, T<:Colorant}
    if wrapped(s)
        error("bond_nonpbc does not support wrapped system. Use bond_pbc instead.")
    end
    bond_radius = bdata[Single].scale[][1][1]
    bond_types = keys(bdata)

    for btype in bond_types
        n = 2 * count(==(btype), bondtypelist)
        resize!(bdata[btype].origin.val   , n)
        resize!(bdata[btype].direction.val, n)
        resize!(bdata[btype].scale.val    , n)
        resize!(bdata[btype].colors.val   , n)
    end

    ibond = Dict{BondTypeVis, Int}(btype => 0 for btype in bond_types)
    for (i, edge) in enumerate(bonds(s))
        atom_id1, atom_id2 = src(edge), dst(edge)
        btype = bondtypelist[i]
        ibond[btype] += 1

        p1, p2 = position(s, atom_id1), position(s, atom_id2)
        center = p1 + 0.5*(p2 - p1)
        bond_length = 0.5*norm(p2 - p1)
        blinv = 1 / bond_length

        bd = bdata[btype]
        n = 2*ibond[btype] - 1
        bd.origin.val[n] = p1
        bd.direction.val[n] = (center - p1) * blinv
        bd.scale.val[n] = Makie.Vec3f(bond_radius, bond_radius, bond_length)
        bd.colors.val[n] = colors[atom_id1]

        bd.origin.val[n+1] = p2
        bd.direction.val[n+1] = (center - p2) * blinv
        bd.scale.val[n+1] = Makie.Vec3f(bond_radius, bond_radius, bond_length)
        bd.colors.val[n+1] = colors[atom_id2]
    end

    return bdata
end

function bond_pbc!(
    bdata::Dict{BondTypeVis, BondData},
    s::AbstractSystem{3, F, SysType},
    colors::AbstractVector{T},
    bondtypelist::Vector{BondTypeVis}
) where {F<:AbstractFloat, SysType<:AbstractSystemType, T<:Colorant}
    if !wrapped(s)
        error("bond_pbc does not support non-wrapped system. Use bond_nonpbc instead.")
    end
    bond_radius = bdata[Single].scale[][1][1]
    bond_types = keys(bdata)

    for btype in bond_types
        n = 5 * count(==(btype), bondtypelist)
        resize!(bdata[btype].origin.val   , n)
        resize!(bdata[btype].direction.val, n)
        resize!(bdata[btype].scale.val    , n)
        resize!(bdata[btype].colors.val   , n)
    end

    origin = box(s).origin # box origin
    axis = box(s).axis
    a, b, c = axis[:, 1], axis[:, 2], axis[:, 3] # box vector
    ibond = Dict{BondTypeVis, Int}(btype => 0 for btype in bond_types)
    for (i, edge) in enumerate(bonds(s))
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
            α = diff[1] == 1 ? 1 : 0
            _t = _cross_point_x(α, p1, bond_vector, a, b, c, origin)
            cross_points[n] = (progress=_t, cross_dim=1)
            n += 1
        end
        if diff[2] != 0
            β = diff[2] == 1 ? 1 : 0
            _t = _cross_point_y(β, p1, bond_vector, a, b, c, origin)
            cross_points[n] = (progress=_t, cross_dim=2)
            n += 1
        end
        if diff[3] != 0
            γ = diff[3] == 1 ? 1 : 0
            _t = _cross_point_z(γ, p1, bond_vector, a, b, c, origin)
            cross_points[n] = (progress=_t, cross_dim=3)
            n += 1
        end
        sort!(cross_points; by=x->x.progress)
        @assert all(x -> -1 <= x.progress <= 1, cross_points)


        # now bond_vector is separated by cross_points [0, ..., 0.5, ..., 1.0]
        color1, color2 = colors[atom_id1], colors[atom_id2]
        btype = bondtypelist[i]
        bd = bdata[btype]
        ibond[btype] += 1
        pbc_shift = SVector{3, F}(0.0, 0.0, 0.0)
        jmax = 1
        for j in 1:length(cross_points)-1
            # bond part begins at t, ends at t_next
            t, t_next = cross_points[j].progress, cross_points[j+1].progress
            # starting point of bond part separated by cross_points
            p  = p1 + t_next*bond_vector + pbc_shift
            bond_vec = (t-t_next)*bond_vector
            bond_len = norm(bond_vec)
            inv = 1 / bond_len

            # add bond part to BondData
            @assert j ≤ 5
            _i = 5*(ibond[btype]-1) + j
            bd.origin.val[_i] = p + bond_len * bond_vec
            bd.direction.val[_i] = inv * bond_vec
            bd.scale.val[_i] = Makie.Vec3f(bond_radius, bond_radius, bond_len)
            bd.colors.val[_i] = t < 0.5 ? color1 : color2

            dim = cross_points[j+1].cross_dim
            if dim != 0
                pbc_shift -= diff[dim] * axis[:, dim]
            end
            jmax = j
        end
        for untouched in jmax+1:5
            _i = 5*(ibond[btype]-1) + untouched
            bd.colors.val[_i] = palette[:transparent]
        end
    end

    return bonds
end

const _cross_point_x, _cross_point_y, _cross_point_z = begin
    @variables p₁ p₂ p₃ # position of bond origin
    @variables v₁ v₂ v₃ # bond vector
    @variables a₁ a₂ a₃ b₁ b₂ b₃ c₁ c₂ c₃ # box vector
    @variables o₁ o₂ o₃ # box origin
    @variables α β γ # coefficients of box vector where αa + βb + γc
    @variables t # progress along bond

    # p + t*v  ~ α*a + β*b + γ*c かつ {α|β|γ} == {1|-1}のときbondが周期境界と交差
    eq₁ = p₁ + t*v₁ ~ α*a₁ + β*b₁ + γ*c₁ + o₁
    eq₂ = p₂ + t*v₂ ~ α*a₂ + β*b₂ + γ*c₂ + o₂
    eq₃ = p₃ + t*v₃ ~ α*a₃ + β*b₃ + γ*c₃ + o₃

    # p, v, a, b, c is given in all cases below
    # if bond crosses yz pbc plane, α == {1|-1}, solve for t, β, γ
    result = Symbolics.solve_for([eq₁, eq₂, eq₃], [t, β, γ]; simplify=true)
    fx = build_function(
        result[1],
        #α, [p₁, p₂, p₃], [v₁, v₂, v₃], [a₁, a₂, a₃], [b₁, b₂, b₃], [c₁, c₂, c₃], [o₁, o₂, o₃];
        α, (p₁, p₂, p₃), (v₁, v₂, v₃), (a₁, a₂, a₃), (b₁, b₂, b₃), (c₁, c₂, c₃), (o₁, o₂, o₃);
        expression=Val{true}
    )

    # if bond crosses xz pbc plane, β == {1|-1}, solve for t, α, γ
    result = Symbolics.solve_for([eq₁, eq₂, eq₃], [t, α, γ]; simplify=true)
    fy = build_function(
        result[1],
        β, [p₁, p₂, p₃], [v₁, v₂, v₃], [a₁, a₂, a₃], [b₁, b₂, b₃], [c₁, c₂, c₃], [o₁, o₂, o₃];
        expression=Val{true}
    )

    # if bond crosses xy pbc plane, γ == {1|-1}, solve for t, α, β
    result = Symbolics.solve_for([eq₁, eq₂, eq₃], [t, α, β]; simplify=true)
    fz = build_function(
        result[1],
        γ, [p₁, p₂, p₃], [v₁, v₂, v₃], [a₁, a₂, a₃], [b₁, b₂, b₃], [c₁, c₂, c₃],[o₁, o₂, o₃];
        expression=Val{true}
    )

    eval(fx), eval(fy), eval(fz)
end
