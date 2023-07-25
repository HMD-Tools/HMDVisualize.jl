const _Bond_Quality = 8
const _Single_Radius = 1.0
const _Double_Radius = 1/3
const _Triple_Radius = 1/5

# origin, extremity, radius1, radius2, segments
const Single_Bond = Makie._mantle(Point3f(zeros(3)), Point3f(0,0,1), _Single_Radius, _Single_Radius, _Bond_Quality)

const Double_Bond = begin
    δ = 0.02
    m1 = Makie._mantle(Point3f(-2/3+δ, 0.0, 0.0), Point3f(-2/3+δ, 0.0, 1.0), _Double_Radius+δ, _Double_Radius+δ, _Bond_Quality)
    m2 = Makie._mantle(Point3f( 2/3-δ, 0.0, 0.0), Point3f( 2/3-δ, 0.0, 1.0), _Double_Radius+δ, _Double_Radius+δ, _Bond_Quality)
    merge([m1, m2])
end

const Triple_Bond = begin
    m1 = Makie._mantle(Point3f(-4/5, 0.0, 1.0), Point3f(-4/5, 0.0, 1.0), _Triple_Radius, _Triple_Radius, _Bond_Quality)
    m2 = Makie._mantle(Point3f( 4/5, 0.0, 1.0), Point3f( 4/5, 0.0, 1.0), _Triple_Radius, _Triple_Radius, _Bond_Quality)
    m3 = Makie._mantle(Point3f( 0.0, 0.0, 0.0), Point3f( 0.0, 0.0, 1.0), _Triple_Radius, _Triple_Radius, _Bond_Quality)
    merge([m1, m2, m3])
end

const Aromatic_Bond = begin
    #m = Makie._mantle(Point3f(0.5, 0.0, 0.0), Point3f(0.5, 0, 1), 0.5, 0.5, _Bond_Quality)
    m = Makie._mantle(Point3f(0.3, 0.0, 0.0), Point3f(0.3, 0.0, 1), 0.65, 0.65, _Bond_Quality)
    n = 2 # number of dots
    interval = [((i-1)/2n, i/2n) for i in 1:2n if isodd(i)]
    dots = [
        #Makie._mantle(Point3f(-3/4, 0.0, start), Point3f(-3/4, 0, finish), 0.25, 0.25, _Bond_Quality)
        Makie._mantle(Point3f(-1.05, 0.0, start), Point3f(-1.05, 0, finish), 0.3, 0.3, _Bond_Quality)
        for (start, finish) in interval
    ]
    merge([m, dots...])
end
