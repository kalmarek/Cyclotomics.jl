using Test
using Statistics
using BenchmarkTools

using Cyclotomics

function zumbroich_perf(ran)
    @testset "zumbroich perf" begin
        @test all(
            Cyclotomics.zumbroich_plain.(ran) .==
            sort!.(collect.(Cyclotomics.zumbroich_basis.(ran))) .==
            Cyclotomics.zumbroich_direct.(ran),
        )

        plain = [(@timed sum(first, Cyclotomics.zumbroich_plain(i)))[2] for i in ran]
        @info "By plain loop: mean and variance over $ran" μ =
            round(mean(plain), sigdigits = 4) σ =
            round(Statistics.var(plain), sigdigits = 4)
        @btime Cyclotomics.zumbroich_plain($(last(ran)))

        compl = [(@timed sum(first, Cyclotomics.zumbroich_basis(i)))[2] for i in ran]
        @info "Constructing the complement: mean and variance over $ran" μ =
            round(mean(compl), sigdigits = 4) σ =
            round(Statistics.var(compl), sigdigits = 4)
        @btime Cyclotomics.zumbroich_basis($(last(ran)))

        direc = [(@timed length(Cyclotomics.zumbroich_direct(i)))[2] for i in ran]
        @info "GAP style: mean and variance over $ran" μ = round(mean(direc), sigdigits = 4) σ =
            round(Statistics.var(direc), sigdigits = 4)
        @btime Cyclotomics.zumbroich_direct($(last(ran)))
    end
end

zumbroich_perf(1:1000)
zumbroich_perf(25000:30001)
