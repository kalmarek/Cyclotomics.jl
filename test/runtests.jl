using Test
using Cyclotomics

@testset "Cyclotomics" begin
    @testset "zumbroich" begin
        @test Cyclotomics._forbidden_residues(45, 5, 1) == BitSet([0])
        @test Cyclotomics._forbidden_residues(45, 3, 2) == BitSet([4, 0, 5])

        @test last(first(Cyclotomics.ForbiddenResidues(9))) isa BitSet
        fb = Cyclotomics.ForbiddenResidues(9)
        @test 0 in fb
        @test 1 in fb
        @test !(2 in fb)

        cs(x) = sort!(collect(x))

        @test cs(Cyclotomics.zumbroich_basis(9)) == [2, 3, 4, 5, 6, 7]

        @test !any(in(fb), Cyclotomics.zumbroich_basis(9))

        @test Cyclotomics.zumbroich_plain(8) ==
              cs(Cyclotomics.zumbroich_basis(8)) ==
              [0, 1, 2, 3]
        @test Cyclotomics.zumbroich_plain(9) ==
              cs(Cyclotomics.zumbroich_basis(9)) ==
              [2, 3, 4, 5, 6, 7]

        @test cs(Cyclotomics.zumbroich_basis(45)) == [
            1,
            2,
            3,
            6,
            7,
            8,
            11,
            12,
            16,
            17,
            19,
            21,
            24,
            26,
            28,
            29,
            33,
            34,
            37,
            38,
            39,
            42,
            43,
            44,
        ]

        @test all(
            Cyclotomics.zumbroich_plain(i) ==
            cs(Cyclotomics.zumbroich_basis(i)) ==
            Cyclotomics.zumbroich_direct(i) for i in 1:5000
        )
    end

    @testset "elementary ops" begin
        @test E(5) isa Number
        @test E(5) isa Cyclotomic{Int}
        @test E(5, 0) isa Cyclotomic{Int}
        @test E(5, 6) isa Cyclotomic{Int}

        E₅ = E(5)
        @test valtype(E₅) == Int
        @test Cyclotomic{Float64}(E₅) isa Cyclotomic{Float64}

        E₅fl = Cyclotomic{Float64}(E₅)
        @test valtype(E₅fl) == Float64

        @test similar(E₅) isa Cyclotomic{Int}
        @test similar(E₅fl) isa Cyclotomic{Float64}

        @test zero(E₅) isa Cyclotomic{Int}
        @test zero(E₅fl) isa Cyclotomic{Float64}

        @test one(E₅) isa Cyclotomic{Int}
        @test one(E₅fl) isa Cyclotomic{Float64}

        @test deepcopy(E₅).coeffs !== E₅.coeffs
        @test deepcopy(E₅).coeffs == E₅.coeffs
    end

    @testset "io strings" begin
        @test sprint(print, E(5)) == " 1*E(5)^1"
        @test sprint(show, 2E(5)) == " 2*ζ₅"
        @test sprint(print, -E(5)) == "-1*E(5)^1"
        @test sprint(show, -E(5)) == "-ζ₅"

        @test sprint(show, -1.0 * E(5)) == "-1.0*ζ₅"

        @test sprint(show, 0.0 * E(4)) == "0.0"
        @test sprint(print, 0.0 * E(4)) == "0.0"
        @test sprint(show, E(1)) == "1"
        @test sprint(print, 1.0 * E(1)) == " 1.0*E(1)^0"

        using Base.Meta
        x = E(5) + 2E(5)^2
        @test sprint(print, x) == " 1*E(5)^1 + 2*E(5)^2"
        @test sprint(show, x) == " ζ₅ + 2*ζ₅²"
        @test eval(Base.Meta.parse(sprint(print, x))) == x

        x = E(5) - 2E(5)^2
        @test sprint(print, x) == " 1*E(5)^1 -2*E(5)^2"
        @test sprint(show, x) == " ζ₅ -2*ζ₅²"
        @test eval(Base.Meta.parse(sprint(print, x))) == x

        x = -2E(5) + 2E(5)^2
        @test sprint(print, x) == "-2*E(5)^1 + 2*E(5)^2"
        @test sprint(show, x) == "-2*ζ₅ + 2*ζ₅²"
        @test eval(Base.Meta.parse(sprint(print, x))) == x
    end

    @testset "indexing and iteration" begin
        x = E(5)
        @test x[0] == 0
        @test x[1] == 1
        @test all(iszero, x[2:5])
        @test x[-1] == 0
        @test x[-4] == 1
        @test x[6] == 1

        @test setindex!(x, 1, 2) == 1
        x[2] = 3
        @test x[0] == 0
        @test x[1] == 1
        @test x[2] == 3
        @test x[-3] == 3
        @test x[7] == 3

        @test collect(Cyclotomics.exps_coeffs(E(6))) == [(1, 1)]
        @test collect(Cyclotomics.exps_coeffs(Cyclotomics.normalform!(E(6)))) ==
              [(4, -1)]
    end

    @testset "aritmetic: +, -, module: *, //" begin
        x = E(5)
        y = E(5, 2)

        @test 2x isa Cyclotomic{Int}
        @test 2.0x isa Cyclotomic{Float64}
        @test x * 2.0 isa Cyclotomic{Float64}
        @test div(x, 2) isa Cyclotomic{Int}
        @test x // 2 isa Cyclotomic{Rational{Int}}
        @test x / 2.0 isa Cyclotomic{Float64}
        @test x / 2 isa Cyclotomic{Float64}

        @test x + 2y isa Cyclotomic{Int}
        xy = x + 2y
        @test xy[0] == 0
        @test xy[1] == 1
        @test xy[2] == 2
        @test all(iszero, xy[3:5])

        @test 2.0xy[0] isa Float64
        @test 2.0xy[1] == 2

        @test (xy-2y)[1] == 1
        @test all(iszero, (xy-2y)[2:5])

        @test 1 + x isa Cyclotomic{Int}
        @test (1+x)[0] == 1
        @test x + 1 isa Cyclotomic{Int}
        @test 2.0 + x isa Cyclotomic{Float64}
        @test 2.0 - x isa Cyclotomic{Float64}
        @test x + 2.0 isa Cyclotomic{Float64}
        @test (x+2.0)[0] == 2.0

        # Bug: normalform! is needed in div
        x = E(4, 0) - E(4, 2)
        @test isone(div(x, 2))
        x = Cyclotomics.embed(5 * E(4, 0), 60)
        Cyclotomics.normalform!(x)
        @test isone(div(x - 2, 3))

        # broadcasting on 1.6 is broken

        @test E(3) .* [1, 2] == [E(3), 2E(3)]
        @test eltype(E(3) .* [1.0, 2.0]) <: Cyclotomic{Float64}
        @test eltype(1 // 1 * E(3) .* [1, 2]) <: Cyclotomic{Rational{Int}}
    end

    @testset "*, powering" begin
        x = E(5)
        y = E(5, 2)

        w = x * (x + 2y)
        @test x * y isa Cyclotomic{Int}
        @test (x*y)[1] == 0
        @test (x*y)[2] == 0
        @test (x*y)[3] == 1
        w = (1 + x) * (1 + y)
        @test w[0] == 1
        @test w[1] == 1
        @test w[2] == 1
        @test w[3] == 1

        @test (x^2)[1] == 0
        @test (x^2)[2] == 1

        @test ((1+x)^2)[0] == 1
        @test ((1+x)^2)[1] == 2
        @test ((1+x)^2)[2] == 1
        @test ((1+x)^2)[3] == 0

        @test ((1+x^3)^2)[0] == 1
        @test ((1+x^3)^2)[1] == 1
        @test ((1+x^3)^2)[3] == 2
    end

    @testset "normal form" begin
        x = E(45) + E(45, 5)
        x.coeffs
        @test deepcopy(x) !== x

        y = Cyclotomics.normalform(x)

        e = E(45)
        w = e + e^2 + e^8 + e^11 + e^17 + e^26 + e^29 + e^38 + e^44

        @test w == y
        @test x.coeffs != y.coeffs
        @test y.coeffs == w.coeffs

        @test deepcopy(x) == deepcopy(y)
        @test x !== deepcopy(x)
        @test x.coeffs === copy(x).coeffs

        @test hash(deepcopy(x)) == hash(deepcopy(y))
        @test length(
            Set([deepcopy(x), deepcopy(y), deepcopy(x), deepcopy(y)]),
        ) == 1

        @test iszero(1 + x - x - 1)

        @test isone(sum(-E(5)^i for i in 1:4))
        @test isone(E(5)^5)
        x = E(5)^5
        @test x == sum(-E(5)^i for i in 1:4)
    end

    @testset "predicates" begin
        @test isreal(1 + E(5) - E(5))
        @test isreal(E(5, 1) + E(5, 4))
        @test isreal(E(5, 2) + E(5, 3))
        @test !isreal(E(5, 1) + E(5, 2))
        @test !isreal(E(5, 1) + E(5, 3))

        @test isreal(abs2(E(5, 1) + E(5, 2)))

        @test 1 == E(5)^5
        @test E(5)^2 + E(5)^3 ≈ (-1 - sqrt(5)) / 2
        @test E(5)^2 + E(5)^3 != (-1 - sqrt(5)) / 2
        @test float(E(5)^2 + E(5)^3) == (-1 - sqrt(5)) / 2
        @test 2.0(E(5)^2 + E(5)^3) ≈ (-1 - sqrt(5))
        @test 2.0(E(5)^2 + E(5)^3) != (-1 - sqrt(5))

        @test isapprox(1e-17E(5), 0.0, atol = 1e-12)
        @test isapprox(0.0, 1e-17E(5), atol = 1e-12)
    end

    @testset "embedding" begin
        let x = E(45)^5 + E(45)^10
            @test conductor(Cyclotomics.reduced_embedding(x)) == 9
            y = Cyclotomics.reduced_embedding(x)
            @test y == E(9)^2 - E(9)^4 - E(9)^7
            Cyclotomics.normalform!(x)
            @test conductor(Cyclotomics.reduced_embedding(x)) == 9
        end

        @test conductor(Cyclotomics.reduced_embedding(E(6))) == 3
        @test conductor(Cyclotomics.reduced_embedding(E(14))) == 7
        @test conductor(Cyclotomics.reduced_embedding(E(1234))) == 617
    end

    @testset "tests against GAP" begin
        @test E(9) == -E(9)^4 - E(9)^7
        @test E(9)^3 == E(3)
        @test E(6) == -E(3)^2
        @test E(12) // 3 == -1 // 3 * E(12)^7

        @test E(45)^4 == -E(45)^19 - E(45)^34
        @test E(45)^13 == -E(45)^28 - E(45)^43
        @test E(45)^14 == -E(45)^29 - E(45)^44
        @test E(45)^22 == -E(45)^7 - E(45)^37

        @test E(5) + E(3) ==
              -E(15)^2 - 2 * E(15)^8 - E(15)^11 - E(15)^13 - E(15)^14
        @test (E(5) + E(5)^4)^2 == -2 * E(5) - E(5)^2 - E(5)^3 - 2 * E(5)^4
        @test E(5) / E(3) == E(15)^13
        @test E(5) * E(3) == E(15)^8
    end

    @testset "conjugation and inverse" begin
        function rand1(α::Cyclotomic, u::AbstractRange, k = 5)
            x = zero(eltype(u)) * α
            for (idx, c) in zip(rand(0:conductor(α), k), rand(u, k))
                x[idx] = c
            end
            return x
        end

        for x in [E(45) + E(45)^2, E(45) + E(45)^2 // 1]
            @test conj(x, 1) == x
            y = prod(
                conj(x, i) for i in 2:conductor(x) if gcd(conductor(x), i) == 1
            )
            @test isreal(x * y)
            @test y ==
                  E(45)^2 + E(45)^3 - E(45)^6 - E(45)^8 + E(45)^11 - E(45)^12 -
                  2 * E(45)^16 +
                  E(45)^17 +
                  E(45)^19 +
                  E(45)^21 - 2 * E(45)^24 - E(45)^26 - E(45)^28 + 2 * E(45)^29 -
                  E(45)^34 + E(45)^37 - 2 * E(45)^42 - E(45)^43 + E(45)^44

            @test Cyclotomics.galois_conj(x, 1) == x
            @test_throws AssertionError Cyclotomics.galois_conj(x, 5)
        end

        for x in [
            E(45)^5 + E(45)^10,
            E(45) - E(45)^5,
            rand1(E(45), -5:5, 3),
            rand1(E(45), -1:1, 5),
        ]
            iszero(x) || @test isone(x * inv(x // big(1)))
        end

        for x in
            [E(45)^5 + big(1) * E(45)^10, (E(45)^5) // 1 + big(1) * E(45)^10]
            y = Cyclotomics.normalform(x)
            @test inv(y) == -E(9)^2 + E(9)^3 - E(9)^4 == inv(x)
            @test inv(y) * x == inv(x) * x == one(x)
        end

        for x in
            [E(45)^5 + big(1) * E(45)^10, (E(45)^5) // 1 + big(1) * E(45)^10]
            y = Cyclotomics.reduced_embedding(x)
            @test inv(y) == -E(9)^2 + E(9)^3 - E(9)^4 == inv(x)
            @test inv(y) * x == inv(x) * x == one(x)
        end

        for x in [0.5 + 0.75 * E(4), 1 // 2 + 3 // 4 * E(4)]
            @test real(x) == 0.5
            @test imag(x) == 0.75
            @test float(real(x)) isa Float64
            @test float(imag(x)) isa Float64
            @test x == real(x) + im * imag(x)
            @test_throws InexactError float(x)
            @test_throws InexactError Rational{Int}(x)
            @test Rational{Int}(x + conj(x)) isa Rational{Int}
            @test Rational{Int}(x + conj(x)) == x + conj(x)
            if valtype(x) <: Rational
                @test Rational(x + conj(x)) isa Rational{Int}
                @test Rational(x + conj(x)) == 1 // 1
            else
                @test_throws MethodError Rational(x + conj(x))
            end
        end

        let x = E(45)^5 + E(45)^10
            @test isreal(x + conj(x))
            @test float(x + conj(x)) isa Float64
            @test !first(Cyclotomics._isreal(x))
            @test_throws InexactError Float64(x)
            @test_throws InexactError Rational(x)
        end

        let x =
                0.9999999999999787 * E(5)^1 +
                0.999999999999992 * E(5)^2 +
                0.9999999999999889 * E(5)^3 +
                0.9999999999999796 * E(5)^4
            z = Cyclotomics.roundcoeffs!(deepcopy(x), digits = 12)
            @test x != -1.0
            @test z == -1.0
            @test inv(z) == -1.0
            @test Cyclotomics.droptol!(inv(x) - inv(z), 1e-12) == 0
        end
    end

    @testset "Conversions" begin
        Cyc = typeof(E(3))

        @test iszero(Cyc(0))
        @test isone(Cyc(1))
        @test isone(Cyc(1.0))
        @test valtype(Cyc(1.0)) == Int
        @test zeros(typeof(E(3)), 2, 2) isa Matrix{<:Cyclotomic}

        v = [E(3)^i for i in 1:3]

        @test (v[1] = 1.0 * E(5)) isa Cyclotomic{Float64}
        @test eltype(v) <: Cyclotomic{Int}
        @test v[1] isa Cyclotomic{Int}
        @test (v[1] = 2.0 * E(5)) isa Cyclotomic{Float64}
        @test v[1] == 2E(5)
        @test_throws InexactError v[1] = 2.5 * E(5)
    end

    @testset "Conversions to julia types" begin
        x = E(3)
        y = x + x^2

        @test Int(y) == -1
        @test y == -1
        @test_throws InexactError Int(x)

        @test float(y) == -1.0
        @test y == -1.0
        @test_throws InexactError float(x)

        @test Float64(y) == -1.0
        # @test y == -1.0
        @test_throws InexactError Float64(x)

        @test float(y // big(1)) isa BigFloat
        @test Float64(y // big(1)) isa Float64

        @test ComplexF64(x) isa ComplexF64
        @test Complex{BigFloat}(x) isa Complex{BigFloat}
        γ = 2 * π / 3
        bγ = 2 * big(π) / 3
        @test Complex{BigFloat}(x) ≈ cos(bγ) + im * sin(bγ)

        @test ComplexF64(x)^2 ≈ ComplexF64(x^2)

        @test complex(x) isa ComplexF64
        @test complex(x) == cos(γ) + im * sin(γ)

        @test complex(big(1) * x) isa Complex{BigFloat}
        @test complex(big(1) * x) ≈ cos(bγ) + im * sin(bγ)

        @test abs(x) ≈ 1.0
        @test abs(big(1) / 3 * x) ≈ big(1) / 3

        @test ComplexF64 == typeof(@inferred ComplexF64(x))
        @test ComplexF64 == typeof(@inferred ComplexF64(y))
        @test ComplexF64 == typeof(@inferred ComplexF64(big(1) * x))
        @test ComplexF64 == typeof(@inferred ComplexF64(big(1) * y))

        @test_logs (
            :error,
            "The cyclotomic is real but it can not be converted to Rational:  1*E(5)^2 + 1*E(5)^3 ≈ -1.618033988749895",
        ) try
            Rational(E(5)^2 + E(5)^3)
        catch ex
            @test ex isa InexactError
        end
    end

    @testset "Complex arithmetic" begin
        @test real(E(4)) == 0
        @test imag(-E(4)) == -1
        @test Cyclotomic(1 + 2im) isa Cyclotomic
        @test Cyclotomic(1.0 - 2.0im) isa Cyclotomic
        @test reim(Cyclotomic(2 - 3im)) == (2, -3)

        @test E(4) * im == -1
        @test E(4) - im == 0

        @test E(9) + im // 3 isa Cyclotomic
    end

    @testset "dense/sparse" begin
        x = E(3)
        @test Cyclotomics.dense(x) isa Cyclotomic{Int,<:DenseVector}
        y = Cyclotomics.dense(x)
        @test Cyclotomics.sparse(y) isa
              Cyclotomic{Int,<:Cyclotomics.SparseVector}

        @test coeffs(x) isa Cyclotomics.SparseVector
        @test coeffs(Cyclotomics.sparse(y)) isa Cyclotomics.SparseVector

        @test coeffs(y) isa Vector
        @test coeffs(Cyclotomics.dense(x)) isa Vector

        @test x == y
    end

    @testset "_enable_intermediate_normalization" begin
        function testmat(p)
            ss = [[i, j] for i in 0:p-1 for j in i+1:p-1]
            return [
                (E(p, i' * reverse(j)) - E(p, i' * j)) // p for i in ss, j in ss
            ]
        end

        M = testmat(6)

        val = Cyclotomics._enable_intermediate_normalization()
        Cyclotomics._enable_intermediate_normalization() = true
        @test isone(inv(M) * M)

        Cyclotomics._enable_intermediate_normalization() = false
        @test_throws OverflowError isone(inv(M) * M)

        Cyclotomics._enable_intermediate_normalization() = val
    end
end
