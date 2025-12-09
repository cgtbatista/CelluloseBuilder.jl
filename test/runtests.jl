using CelluloseBuilder
using Test
using TestItemRunner
using TestItems

_stackingfaults_path = joinpath(@__DIR__, "faults/")
_vmd_test_path = joinpath(@__DIR__, "thirdparties/")

@testset "CelluloseBuilder.jl     " begin
    # Write your tests here.
    @test true
end

@testset "Third-party verification" begin
    tclscript = joinpath(_vmd_test_path, "hello.tcl")
    false_exec = tempname()
    Base.open(false_exec, "w") do file
        println(file, "#!/bin/bash")
        println(file, "echo 'False executable!'")
    end
    false_file = tempname() * ".tcl"
    @testset "internal psfgen executable" begin
        ## Is PSFGEN executable?
        result = execPSFGEN(tclscript)
        @test occursin("Hiho, world! VMD/PSFGEN was recognized!", result)
        @test_throws ErrorException execPSFGEN(CelluloseBuilder.psfgen, false_file)
    end
    @testset "Checking executability" begin
        ## When tere is and is not an executable
        @test CelluloseBuilder.hasExecutable(CelluloseBuilder.psfgen)
        @test !(CelluloseBuilder.hasExecutable(false_exec))
    end
    @testset "Skipping false executables" begin
        ## False executable
        @test_throws ErrorException execVMD(false_exec, tclscript)
        @test_throws ErrorException execPSFGEN(false_exec, tclscript)
    end
    @testset "external executables (not so necessary yet)" begin
        ## Is VMD executable?
        if CelluloseBuilder.hasExecutable("vmd")
            result = execVMD("vmd", tclscript)
            @test occursin("Hiho, world! VMD/PSFGEN was recognized!", result) && occursin("Info) VMD for", result)
            @test_throws ErrorException execVMD("vmd", false_file)
        else
            @warn "VMD executable not found. Don't worry, it's not so necessary to default executing yet."
        end
        ## Is PSFGEN executable?
        if CelluloseBuilder.hasExecutable("psfgen") 
            result = execPSFGEN("psfgen", tclscript)
            @test occursin("Hiho, world! VMD/PSFGEN was recognized!", result) && occursin("PSFGEN ", result)
            @test_throws ErrorException execPSFGEN(CelluloseBuilder.psfgen, false_file)
        else
            @warn "PSFGEN executable not found. Don't worry, it's not so necessary to default executing yet."
        end
    end
    Base.rm(false_exec, force=true)
end