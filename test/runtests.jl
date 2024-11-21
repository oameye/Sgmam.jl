using Sgmam
using Test

@testset "Code quality" begin
  using ExplicitImports, Aqua
  ignore_deps = [:Test]

  @test check_no_stale_explicit_imports(Sgmam) == nothing
  @test check_all_explicit_imports_via_owners(Sgmam) == nothing
  Aqua.test_ambiguities(Sgmam)
  Aqua.test_all(
    Sgmam;
    deps_compat = (
        ignore=ignore_deps,
        check_extras=(ignore=ignore_deps,),
        check_weakdeps=(ignore=ignore_deps,)
        ),
    ambiguities = false,
  )
end

@testset "Code linting" begin
  using JET
  JET.test_package(Sgmam; target_defined_modules=true)
end
