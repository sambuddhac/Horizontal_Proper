model = Model(GLPK.Optimizer)
@variable model begin
    x1<=3
    x2<=4
end
@objective (model, Min, X1+X2)
@show(value.(x1),value.(x2))
