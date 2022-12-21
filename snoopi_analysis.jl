using SnoopCompile
using Dynare
#@dynare "test/models/example1/example1.mod"
inf = @snoopi_deep dynare("test/models/example3report/example3report.mod")
@show inf

#using ProfileView 
#ProfileView.view(flamegraph(inf))

#using PProf
#pprof(flamegraph(inf))