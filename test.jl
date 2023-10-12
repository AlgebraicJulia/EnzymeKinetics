
function split(cat,prod1,prod2,on::T) where T
  in = Symbol(first(cat), sep_sym,first(prod1), sep_sym,first(prod2))
  Open(LabelledReactionNet{T,Number}(unique((in=>0, cat,prod1,prod2)), ((Symbol(:split_,in),on),in=>(first(cat),first(prod1),first(prod2)))));
end;

# Beginning of a multisplit motif. Not sure if needed.
function multisplit(enzs,prods,on::T) where T
  in = Symbol(first(cat), sep_sym,first(prod1), sep_sym,first(prod2))
  Open(LabelledReactionNet{T,Number}(unique((in=>0, cat,prod1,prod2)), ((Symbol(:split_,in),on),in=>(first(cat),first(prod1),first(prod2)))));
end;

function split(cat,prod1,prod2)
  in = Symbol(cat, sep_sym,prod1, sep_sym,prod2)
  Open(LabelledPetriNet(unique((in, cat,prod1,prod2)), (Symbol(:split_,in),in=>(cat,prod1,prod2))));
end;

enzXsubYZ = @relation (X, Xinact, Xdeg, YZ, Y, Z) where (X, Xinact, Xdeg, YZ, XYZ, Y, Z) begin
  bindXYZ(X, YZ, XYZ)
  splitXYZ(XYZ, X, Y, Z)
end

# TODO:
# 1. Currently assumes the structures/concentrations are initialize for the end products,
#    rather the substrate that is split. Ie, YZ is set to zero.
#    Change to have YZ first. Perhaps also determine Y and Z from YZ
# 2. Need to determine appropriate bundling of legs. These are original vects.
function enz_sub2(rxns, cat1, sub1, sub2)
  catsym = first(cat1)
  subsym1 = first(sub1)
  subsym2 = first(sub2)
  subsym = Symbol(subsym1, sep_sym, subsym2)
  catsub = Symbol(catsym, sep_sym, subsym)
  obtype = valtype(rates(apex(first(last(first(rxns))))))
  out = oapply(enzXsubYZ, Dict([:bindXYZ, :splitXYZ] .=> rxns[catsub]), Dict(
    :X=>ob(obtype, cat1),
    :Xinact=>ob(obtype, Symbol(catsym,:_inact)=>0),
    :Xdeg=>ob(obtype, Symbol(catsym,:_deg)=>0),
    :YZ=>ob(obtype, subsym=>0),
    :XYZ=>ob(obtype, catsub=>0),
    :Y=>ob(obtype, sub1),
    :Z=>ob(obtype, sub2)))
  bundle_legs(out, [[1,2,3], [4,5]])
end

function enz_sub2(cat1::Symbol, sub1::Symbol, sub2::Symbol)
  catsym = first(cat1)
  subsym1 = first(sub1)
  subsym2 = first(sub2)
  subsym = Symbol(subsym1, sep_sym, subsym2)
  catsub = Symbol(catsym, sep_sym, subsym)
  out = oapply(enzXsubYZ, Dict(:bindXYZ=>bindunbind(cat1, subsym), :splitXYZ=>split(cat1, sub1, sub2)), Dict(
    :X=>ob(cat1),
    :Xinact=>ob(Symbol(catsym,:_inact)),
    :Xdeg=>ob(Symbol(catsym,:_deg)),
    :YZ=>ob(subsym),
    :XYZ=>ob(catsub),
    :Y=>ob(sub1),
    :Z=>ob(sub2)))
  bundle_legs(out, [[1,2,3], [4,5]])
end

# Just old copy at moment. Not sure if needed. Could be used to set up multisplitting
function enzyme_generators(enzymes::Array{Symbol}, substrates::Array{Symbol})
  gens = Dict{Symbol, Any}()
  for e1 in enzymes
    for e2 in enzymes
      if e1 == e2
        gens[Symbol(:cat, e1)] = enz(e1)
      else
        gens[Symbol(:cat, e1, :cat, e2)] = enz_enz(e1, e2)
      end
    end
    for s in substrates
      gens[Symbol(:cat, e1, :sub, s)] = enz_sub(e1, s)
    end
  end
  gens
end
