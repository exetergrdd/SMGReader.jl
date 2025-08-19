# src/test_buffer.jl


@testset "VectorBuffer Construction, length, resize" begin
    buff = VectorBuffer{Int}(10, 5)
    @test buff.len == 5
    @test length(buff.data) == 10

    @test length(buff) == 5
    @test size(buff) == (5,)
    
    @test axes(buff) == (Base.OneTo(5),)
    @test eachindex(buff) == Base.OneTo(5)

    ### shrinking
    setlength!(buff, 3)
    @test length(buff) == 3
    @test length(buff.data) == 10

    ### growing beneath length
    setlength!(buff, 8)
    @test length(buff) == 8
    @test length(buff.data) == 10

     ### growing above length
    setlength!(buff, 20)
    @test length(buff) == 20
    @test length(buff.data) >= 20
end


@testset "VectorBuffer Get/Set Index + iteration + broadcasting" begin
    buff = VectorBuffer{Int}(10, 5)
    buff[1] = 42
    buff[5] = 99
    @test buff[1] == 42
    @test buff[5] == 99
    @test_throws BoundsError buff[0]
    @test_throws BoundsError buff[6]

    fill!(buff, 12)
    @test all(==(12), buff)
    buff .= 2
    @test all(==(2), buff)

    buff .= 1:length(buff)
    inds = [1, 0, 1, 0, 1] .== true
    @test buff[inds] == [1, 3, 5]

    @test buff[2:4] == 2:4
end
