
using DrWatson, Test
@quickactivate :ImmuneBoostingHealthcare

# Here you include files using `srcdir`
# include(srcdir("file.jl"))

# Run test suite
println("Starting tests")
ti = time()

@testset "Model tests" begin
    @testset "Test functions with numbers transitioning" begin
        @testset "Get time for empty model is 0" begin 
            a = HCWSEIRRRVOutput() 
            @test gettime(a) == 0 
        end
        @testset "Empty model has population 0" begin
            a = HCWSEIRRRVOutput()  
            @test get_n(a) == 0
        end   
        @testset "Two empty model outputs equal each other" begin
            a = HCWSEIRRRVOutput()
            b = HCWSEIRRRVOutput() 
            @test a == b
        end
        @testset "Initial time for model with population 10 is 0" begin 
            a = HCWSEIRRRVOutput(10) 
            @test gettime(a) == 0 
        end
        @testset "Model with 10 susceptibles has a population 10" begin
            a = HCWSEIRRRVOutput(10)  
            @test get_n(a) == 10
        end  
        @testset "Model with population 10 differs from model with population 0" begin
            a = HCWSEIRRRVOutput()
            b = HCWSEIRRRVOutput(10) 
            @test a != b
        end
        @testset "Can count 10 susceptibles" begin
            a = HCWSEIRRRVOutput(10)  
            @test get_S(a) == 10
        end  
        @testset "Expose 1 does not change population size" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            @test get_n(a) == 10
        end
        @testset "Expose 1 reduces number susceptible" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            @test get_S(a) == 9
        end
        @testset "Expose 1 leaves 1 exposed" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            @test get_E(a) == 1
        end
        @testset "Expose 1 twice leaves 2 exposed" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            expose!(a, 1) 
            @test get_E(a) == 2
        end
        @testset "Expose more than population generates warning" begin
            a = HCWSEIRRRVOutput(10)  
            @test_logs (:warn, "Negative value in compartment S, -2, at time 0") expose!(a, 12) 
        end
        @testset "No warning if not generating negative value" begin
            a = HCWSEIRRRVOutput(10)  
            @test_nowarn expose!(a, 1)  
        end
        @testset "Advance time increases time by 1" begin
            a = HCWSEIRRRVOutput()
            advancetime!(a)
            @test gettime(a) == 1 
        end  
        @testset "Advance time twice increases time by 2" begin
            a = HCWSEIRRRVOutput()
            advancetime!(a)
            advancetime!(a)
            @test gettime(a) == 2
        end  
        @testset "Expose more than population generates warning after time advanced" begin
            a = HCWSEIRRRVOutput(10)  
            advancetime!(a)
            @test_logs (:warn, "Negative value in compartment S, -2, at time 1") expose!(a, 12) 
        end
        @testset "Progression from exposed to infectious does not change population" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            progress!(a, 1) 
            @test get_n(a) == 10 
        end
        @testset "Progression from exposed to infectious increases number infectious" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            progress!(a, 1) 
            @test get_I(a) == 1 
        end
        @testset "A second progression does not change population" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 2) 
            progress!(a, 1) 
            progress!(a, 1) 
            @test get_n(a) == 10 
        end
        @testset "Two progressions leave 2 infectious" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 2) 
            progress!(a, 1) 
            progress!(a, 1) 
            @test get_I(a) == 2
        end
        @testset "Progression without exposure gives an warning" begin
            a = HCWSEIRRRVOutput(10)  
            @test_logs (:warn, "Negative value in compartment E, -1, at time 0") progress!(a, 1) 
        end
        @testset "Progression following exposure does not give an warning" begin
            a = HCWSEIRRRVOutput(10) 
            expose!(a, 2)  
            @test_nowarn progress!(a, 1) 
        end
        @testset "Subsequent progression without exposure gives progressive warning" begin
            a = HCWSEIRRRVOutput(10)  
            @test_logs (:warn, "Negative value in compartment E, -1, at time 0") progress!(a, 1) 
            @test_logs (:warn, "Negative value in compartment E, -2, at time 0") progress!(a, 1) 
        end
        @testset "Advancing time modifies progression warning" begin
            a = HCWSEIRRRVOutput(10)  
            @test_logs (:warn, "Negative value in compartment E, -1, at time 0") progress!(a, 1) 
            advancetime!(a)
            @test_logs (:warn, "Negative value in compartment E, -2, at time 1") progress!(a, 1) 
        end
        @testset "Diagnosis does not change population size" begin 
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            progress!(a, 1) 
            diagnose!(a, 1) 
            @test get_n(a) == 10
        end
        @testset "Diagnosis reduces number undiagnosed" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            progress!(a, 1) 
            diagnose!(a, 1) 
            @test get_I(a) == 0
        end
        @testset "Diagnosis leaves one newly diagnosed" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            progress!(a, 1) 
            diagnose!(a, 1) 
            @test countnewlydiagnosed(a) == 1
        end
        @testset "Two diagnoses leaves two newly diagnosed" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 2) 
            progress!(a, 2) 
            diagnose!(a, 1) 
            diagnose!(a, 1) 
            @test countnewlydiagnosed(a) == 2
        end
        @testset "Diagnosing more than number undiagnosed gives warning" begin
            a = HCWSEIRRRVOutput(10)  
            @test_logs (:warn, "Negative value in compartment I, -1, at time 0") diagnose!(a, 1) 
        end
        @testset "No warning when diagnosing no more than number undiagnosed" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 2) 
            progress!(a, 2) 
            @test_nowarn diagnose!(a, 1) 
        end
        @testset "Diagnosing more than undiagnosed twice gives progressive warning" begin
            a = HCWSEIRRRVOutput(10)  
            @test_logs (:warn, "Negative value in compartment I, -1, at time 0") diagnose!(a, 1) 
            @test_logs (:warn, "Negative value in compartment I, -2, at time 0") diagnose!(a, 1) 
        end
        @testset "Diagnosis warning responds appropriately to advancing time " begin
            a = HCWSEIRRRVOutput(10)  
            @test_logs (:warn, "Negative value in compartment I, -1, at time 0") diagnose!(a, 1) 
            advancetime!(a)
            @test_logs (:warn, "Negative value in compartment I, -2, at time 1") diagnose!(a, 1) 
        end
        @testset "One diagnosis gives total of 1 diagnosed" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            progress!(a, 1) 
            diagnose!(a, 1) 
            @test gettotaldiagnosed(a) == 1
        end
        @testset "Two diagnoses gives total of 2 diagnosed" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 2) 
            progress!(a, 2) 
            diagnose!(a, 2) 
            @test gettotaldiagnosed(a) == 2
        end
        @testset "Advancing time reduces number newly diagnosed" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            progress!(a, 1) 
            diagnose!(a, 1) 
            advancetime!(a)
            @test countnewlydiagnosed(a) == 0
        end
        @testset "Advancing time does not change population" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            progress!(a, 1) 
            diagnose!(a, 1) 
            advancetime!(a)
            @test get_n(a) == 10
        end
        @testset "After advancing time only newly diagnosed are counted as such" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 4) 
            progress!(a, 4) 
            diagnose!(a, 2) 
            advancetime!(a)
            diagnose!(a, 2) 
            @test countnewlydiagnosed(a) == 2
        end
        @testset "Advancing time twice with diagnoses in between does not change population" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 4) 
            progress!(a, 4) 
            diagnose!(a, 2) 
            advancetime!(a)
            diagnose!(a, 2) 
            advancetime!(a)
            @test get_n(a) == 10
        end
        @testset "Advancing time twice with diagnoses in between gives sum of diagnoses as total" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 4) 
            progress!(a, 4) 
            diagnose!(a, 2) 
            advancetime!(a)
            diagnose!(a, 2) 
            advancetime!(a)
            @test gettotaldiagnosed(a) == 4
        end
        @testset "Advancing time nine times does not change population" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            progress!(a, 1) 
            diagnose!(a, 1) 
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            @test get_n(a) == 10
        end
        @testset "Advancing time nine times does not reduce number diagnosed" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            progress!(a, 1) 
            diagnose!(a, 1) 
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            @test gettotaldiagnosed(a) == 1
        end
        @testset "Advancing time ten times does not change population" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            progress!(a, 1) 
            diagnose!(a, 1) 
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            @test get_n(a) == 10
        end
        @testset "Advancing time nine times does reduce number diagnosed" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            progress!(a, 1) 
            diagnose!(a, 1) 
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            @test gettotaldiagnosed(a) == 0
        end
        @testset "Advancing time nine times makes 1 immune" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            progress!(a, 1) 
            diagnose!(a, 1) 
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            @test gettotalimmune(a) == 1
        end
        @testset "Waning from R1 does not change population size" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            progress!(a, 1) 
            diagnose!(a, 1) 
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            wane1!(a, 1)
            @test get_n(a) == 10
        end
        @testset "Waning from R1 does not change number immune" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            progress!(a, 1) 
            diagnose!(a, 1) 
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            wane1!(a, 1)
            @test gettotalimmune(a) == 1
        end
        @testset "Waning more than are in R1 generates warning" begin
            a = HCWSEIRRRVOutput(10)  
            @test_logs (:warn, "Negative value in compartment R1, -1, at time 0") wane1!(a, 1)
        end
        @testset "Waning from R2 does not change population size" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            progress!(a, 1) 
            diagnose!(a, 1) 
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            wane1!(a, 1)
            wane2!(a, 1)
            @test get_n(a) == 10
        end
        @testset "Waning from R2 does not change total number immune" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 2) 
            progress!(a, 2) 
            diagnose!(a, 2) 
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            wane1!(a, 2)
            wane2!(a, 2)
            @test gettotalimmune(a) == 2
        end
        @testset "Waning from R3 does not change population size" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            progress!(a, 1) 
            diagnose!(a, 1) 
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            wane1!(a, 1)
            wane2!(a, 1)
            wane3!(a, 1)
            @test get_n(a) == 10
        end
        @testset "Waning from R3 does change total number immune" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 4) 
            progress!(a, 4) 
            diagnose!(a, 4) 
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            wane1!(a, 3)
            wane2!(a, 2)
            wane3!(a, 1)
            @test gettotalimmune(a) == 3
        end
        @testset "Vaccinating susceptibles does not change population" begin
            a = HCWSEIRRRVOutput(10)  
            vaccinateS!(a, 2)
            @test get_n(a) == 10
        end
        @testset "Vaccinating susceptibles increases number vaccine pending" begin
            a = HCWSEIRRRVOutput(10)  
            vaccinateS!(a, 2)
            @test get_v(a) == 2
        end
        @testset "Movement from vaccine pending does not change population" begin
            a = HCWSEIRRRVOutput(10)  
            vaccinateS!(a, 2)
            becomeimmunefromvaccine!(a, 2)
            @test get_n(a) == 10
        end
        @testset "Movement from vaccine pending does change total number immune" begin
            a = HCWSEIRRRVOutput(10)  
            vaccinateS!(a, 2)
            becomeimmunefromvaccine!(a, 2)
            @test gettotalimmune(a) == 2
        end
        @testset "Exposure after vaccine does not change population" begin
            a = HCWSEIRRRVOutput(10)  
            vaccinateS!(a, 2)
            exposevaccinated!(a, 2)
            @test get_n(a) == 10
        end
        @testset "Exposure after vaccine changes number exposed" begin
            a = HCWSEIRRRVOutput(10)  
            vaccinateS!(a, 2)
            exposevaccinated!(a, 2)
            @test get_E(a) == 2
        end
        @testset "Vaccinating from R2 does not change population" begin
            a = HCWSEIRRRVOutput(10)
            expose!(a, 4) 
            progress!(a, 4) 
            diagnose!(a, 4) 
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            wane1!(a, 3)
            vaccinateR2!(a, 2)
            @test get_n(a) == 10
        end
        @testset "Vaccinating from R2 does not change number immune" begin
            a = HCWSEIRRRVOutput(10)
            expose!(a, 4) 
            progress!(a, 4) 
            diagnose!(a, 4) 
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            wane1!(a, 3)
            vaccinateR2!(a, 2)
            @test gettotalimmune(a) == 4
        end
        @testset "Vaccinating from R3 does not change population" begin
            a = HCWSEIRRRVOutput(10)
            expose!(a, 4) 
            progress!(a, 4) 
            diagnose!(a, 4) 
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            wane1!(a, 3)
            wane2!(a, 3)
            vaccinateR3!(a, 2)
            @test get_n(a) == 10
        end
        @testset "Vaccinating from R3 does not change number immune" begin
            a = HCWSEIRRRVOutput(10)
            expose!(a, 4) 
            progress!(a, 4) 
            diagnose!(a, 4) 
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            wane1!(a, 3)
            wane2!(a, 3)
            vaccinateR3!(a, 2)
            @test gettotalimmune(a) == 4
        end
        @testset "Recovery does not change population size" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            progress!(a, 1) 
            recover!(a, 1) 
            @test get_n(a) == 10
        end
        @testset "Recovery increases number immune" begin
            a = HCWSEIRRRVOutput(10)  
            expose!(a, 1) 
            progress!(a, 1) 
            recover!(a, 1) 
            @test gettotalimmune(a) == 1
        end
        @testset "Boosting from R2 does not change population" begin
            a = HCWSEIRRRVOutput(10)
            expose!(a, 4) 
            progress!(a, 4) 
            diagnose!(a, 4) 
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            wane1!(a, 3)
            boostR2!(a, 2)
            @test get_n(a) == 10
        end
        @testset "Boosting from R2 does not change number immune" begin
            a = HCWSEIRRRVOutput(10)
            expose!(a, 4) 
            progress!(a, 4) 
            diagnose!(a, 4) 
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            wane1!(a, 3)
            boostR2!(a, 2)
            @test gettotalimmune(a) == 4
        end
        @testset "Boosting from R3 does not change population" begin
            a = HCWSEIRRRVOutput(10)
            expose!(a, 4) 
            progress!(a, 4) 
            diagnose!(a, 4) 
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            wane1!(a, 3)
            wane2!(a, 3)
            boostR3!(a, 2)
            @test get_n(a) == 10
        end
        @testset "Boosting from R3 does not change number immune" begin
            a = HCWSEIRRRVOutput(10)
            expose!(a, 4) 
            progress!(a, 4) 
            diagnose!(a, 4) 
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            advancetime!(a)
            wane1!(a, 3)
            wane2!(a, 3)
            vaccinateR3!(a, 2)
            @test gettotalimmune(a) == 4
        end
    end
    @testset "Test functions with parameters to guide transitions" begin
        @testset "When sigma == 0.5, half of exposed progress" begin
            p = HCWSEIRRRParameters(; sigma = 0.5)
            a = HCWSEIRRRVOutput(24.0)
            expose!(a, 12) 
            progress!(a, 6) 
            b = HCWSEIRRRVOutput(24.0)
            expose!(b, 12) 
            progress!(b, p) 
            @test a == b 
        end
        @testset "When sigma == 0.5, half of exposed progress each time" begin
            p = HCWSEIRRRParameters(; sigma = 0.5)
            a = HCWSEIRRRVOutput(24.0)
            expose!(a, 12) 
            progress!(a, 6) 
            progress!(a, 3) 
            b = HCWSEIRRRVOutput(24.0)
            expose!(b, 12) 
            progress!(b, p) 
            progress!(b, p) 
            @test a == b 
        end
        @testset "Get an error if sigma < 0" begin
            p = HCWSEIRRRParameters(; sigma=-1.0)
            a = HCWSEIRRRVOutput(24.0)
            expose!(a, 12) 
            @test_throws DomainError progress!(a, p)
        end
        @testset "Get an error if sigma >1" begin
            p = HCWSEIRRRParameters(; sigma=1.5)
            a = HCWSEIRRRVOutput(24.0)
            expose!(a, 12) 
            @test_throws DomainError progress!(a, p)
        end
        @testset "When theta == 0.1, 1/10 of undiagnosed get diagnosed" begin
            p = HCWSEIRRRParameters(; theta=0.1)
            a = HCWSEIRRRVOutput(24)
            expose!(a, 12) 
            progress!(a, 10) 
            diagnose!(a, p)
            @test gettotaldiagnosed(a) == 1 
        end
        @testset "When theta == 0.1, 1/10 of undiagnosed get diagnosed each time" begin
            p = HCWSEIRRRParameters(; theta=0.1)
            a = HCWSEIRRRVOutput(24.0)
            expose!(a, 12.0) 
            progress!(a, 10.0) 
            diagnose!(a, p)
            diagnose!(a, p)
            @test gettotaldiagnosed(a) == 1.9 
        end
        @testset "Get an error if theta < 0" begin
            p = HCWSEIRRRParameters(; theta=-0.1)
            a = HCWSEIRRRVOutput(24.0)
            expose!(a, 12.0) 
            progress!(a, 10.0) 
            @test_throws DomainError diagnose!(a, p)
        end
        @testset "Get an error if theta >1" begin
            p = HCWSEIRRRParameters(; theta=2.0)
            a = HCWSEIRRRVOutput(24.0)
            expose!(a, 12.0) 
            progress!(a, 10.0) 
            @test_throws DomainError diagnose!(a, p)
        end
        @testset "Function _parametermin finds smallest parameter in struct" begin
            p = HCWSEIRRRParameters(; theta=0.1)
            @test ImmuneBoostingHealthcare._parametermin(p) == 0 
        end
        @testset "Function _parametermin finds smallest parameter in struct" begin
            p = HCWSEIRRRParameters(; theta=-0.1)
            @test ImmuneBoostingHealthcare._parametermin(p) == -0.1 
        end
        @testset "When gamma == 0.2, 1/5 of undiagnosed recover" begin
            p = HCWSEIRRRParameters(; gamma = 0.2)
            a = HCWSEIRRRVOutput(24.0)
            expose!(a, 12.0) 
            progress!(a, 10.0) 
            recover!(a, p)
            @test gettotalimmune(a) == 2 
        end
        @testset "If all are diagnosed, none go straight to recovery" begin
            p = HCWSEIRRRParameters(; theta=1.0, gamma = 0.2)
            a = HCWSEIRRRVOutput(24.0)
            expose!(a, 12.0) 
            progress!(a, 10.0) 
            recover!(a, p)
            @test gettotalimmune(a) == 0 
        end
        @testset "When omega = 0.01 then 3 / 100 wane from R1" begin
            p = HCWSEIRRRParameters(; omega=0.01)
            a = HCWSEIRRRVOutput(100.0)
            expose!(a, 100.0) 
            progress!(a, 100.0) 
            recover!(a, 100.0)
            wane1!(a, p)
            @test get_R1(a) == 97 
            @test get_R2(a) == 3 
        end
        @testset "Get an error if omega > 1/3" begin
            p = HCWSEIRRRParameters(; omega=0.4)
            a = HCWSEIRRRVOutput(100.0)
            @test_throws DomainError wane1!(a, p)
        end
        @testset "When omega = 0.01 then 3 / 100 wane from R2" begin
            p = HCWSEIRRRParameters(; omega=0.01)
            a = HCWSEIRRRVOutput(100.0)
            expose!(a, 100.0) 
            progress!(a, 100.0) 
            recover!(a, 100.0)
            wane1!(a, p)
            wane2!(a, p)
            @test get_R2(a) == 2.91 
            @test get_R3(a) == 0.09 
        end
        @testset "When omega = 0.01 then 3 / 100 wane from R3" begin
            p = HCWSEIRRRParameters(; omega=0.01)
            a = HCWSEIRRRVOutput(10_000.0)
            expose!(a, 10_000) 
            progress!(a, 10_000) 
            recover!(a, 10_000 / 3)
            wane1!(a, p)
            wane2!(a, p)
            wane3!(a, p)
            @test get_R3(a) == 2.91 
            @test get_S(a) == 0.09 
        end
        @testset "When nu is a number, vaccinating susceptibles increases vaccine pending" begin
            p = HCWSEIRRRParameters(; nu=0.01)
            a = HCWSEIRRRVOutput(200.0)            
            vaccinateS!(a, p)
            @test get_v(a) == 2
        end
        @testset "When nu is a number, get an error when it is greater than 1" begin
            p = HCWSEIRRRParameters(; nu=1.5)
            a = HCWSEIRRRVOutput(24.0)
            @test_throws DomainError vaccinateS!(a, p)
        end
        @testset "When nu is a function, do not get an error if function returns values in range [0, 1]" begin
            p = HCWSEIRRRParameters(; nu = t -> 0.01)
            a = HCWSEIRRRVOutput(24.0)
            vaccinateS!(a, p, 1)
            @test get_v(a) == 0.24
        end
        @testset "When nu is a function, do get an error if function returns negative vaalue" begin
            p = HCWSEIRRRParameters(; nu = t -> (10 - t) / 10)
            a = HCWSEIRRRVOutput(24.0)
            @test_throws DomainError vaccinateS!(a, p, 11)
        end
    end
end


#=
@testset "Test two empty model outputs equal each other" begin
    m = HCWSEIRRROutput()
    n = HCWSEIRRROutput() 
    @test m == n
end
=#
ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
