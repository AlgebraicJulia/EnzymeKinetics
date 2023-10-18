#=
function split1(cat,prod1,prod2,on::T) where T
  in = Symbol(first(cat), sep_sym,first(prod1), sep_sym,first(prod2))
  Open(LabelledReactionNet{T,Number}(unique((in=>0, cat,prod1,prod2)), ((Symbol(:split_,in),on),in=>(first(cat),first(prod1),first(prod2)))));
end;

function split1(cat,prod1,prod2)
  in = Symbol(cat, sep_sym,prod1, sep_sym,prod2)
  Open(LabelledPetriNet(unique((in, cat,prod1,prod2)), (Symbol(:split_,in),in=>(cat,prod1,prod2))));
end;

# Beginning of a multisplit motif. Not sure if needed.
function multisplit(enzs,prods,on::T) where T
  in = Symbol(first(cat), sep_sym,first(prod1), sep_sym,first(prod2))
  Open(LabelledReactionNet{T,Number}(unique((in=>0, cat,prod1,prod2)), ((Symbol(:split_,in),on),in=>(first(cat),first(prod1),first(prod2)))));
end;


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
=#

function split2(cat,sub,site::Int,on::T) where T
  catsym = first(cat)
  subsym = first(sub)
  in = Symbol(catsym,site,sep_sym,subsym)
  prod1 = Symbol(String(subsym)[1:site])
  prod2 = Symbol(String(subsym)[site+1:end])
  Open(LabelledReactionNet{T,Number}(unique((in=>0, cat,prod1=>0,prod2=>0)), ((Symbol(:split_,in),on),in=>(first(cat),prod1,prod2))));
end;

function split2(cat,sub,site::Int)
  in = Symbol(cat,site,sep_sym,sub)
  prod1 = Symbol(String(sub)[1:site])
  prod2 = Symbol(String(sub)[site+1:end])
  Open(LabelledPetriNet(unique((in, cat,prod1,prod2)), (Symbol(:split_,in),in=>(cat,prod1,prod2))));
end;

function bindunbind_multisite(in1, in2, site::Int, on::T, off::T) where T
  out = Symbol(first(in1),site, sep_sym,first(in2))
  Open(LabelledReactionNet{T,Number}(unique((in1, in2, out=>0)), ((Symbol(:bind,site,:_,first(in1), sep_sym,first(in2)),on),(first(in1),first(in2))=>out),
    ((Symbol(:unbind_,out),off),out=>(first(in1),first(in2)))))
end;

function bindunbind_multisite(in1, in2, site::Int)
  out = Symbol(in1,site, sep_sym,in2)
  Open(LabelledPetriNet(unique((in1, in2,out)), (Symbol(:bind,site,:_,in1, sep_sym,in2),(in1,in2)=>out),
                                                (Symbol(:unbind_,out),out=>(in1,in2))))
end;

enzXsubYZ = @relation (X, Xinact, Xdeg, YZ, Y, Z) where (X, Xinact, Xdeg, YZ, XYZ, Y, Z) begin
  bindXYZ(X, YZ, XYZ)
  splitXYZ(XYZ, X, Y, Z)
end


function enz_sub_split(rxns, cat1, sub, site)
  catsym = first(cat1)
  subsym = first(sub)

  catsub = Symbol(catsym, site, sep_sym, subsym)
  frag1sym = Symbol(String(subsym)[1:site])
  frag2sym = Symbol(String(subsym)[site+1:end])
    
  obtype = valtype(rates(apex(first(last(first(rxns))))))
  out = oapply(enzXsubYZ, Dict([:bindXYZ, :splitXYZ] .=> rxns[catsub]), Dict(
    :X=>ob(obtype, cat1),
    :Xinact=>ob(obtype, Symbol(catsym,:_inact)=>0),
    :Xdeg=>ob(obtype, Symbol(catsym,:_deg)=>0),
    :YZ=>ob(obtype, sub),
    :XYZ=>ob(obtype, catsub=>0),
    :Y=>ob(obtype, frag1sym=>0),
    :Z=>ob(obtype, frag2sym=>0)))
  bundle_legs(out, [[1,2,3], [4,5]])
end

function enz_sub_split(cat1::Symbol, sub::Symbol, site)
  catsym = cat1
  subsym = sub

  catsub = Symbol(catsym, site, sep_sym, subsym)
  frag1sym = Symbol(String(subsym)[1:site])
  frag2sym = Symbol(String(subsym)[site+1:end])

  out = oapply(enzXsubYZ, Dict(:bindXYZ=>bindunbind_multisite(catsym,subsym,site), :splitXYZ=>split2(catsym,subsym,site)), Dict(
    :X=>ob(cat1),
    :Xinact=>ob(Symbol(catsym,:_inact)),
    :Xdeg=>ob(Symbol(catsym,:_deg)),
    :YZ=>ob(sub),
    :XYZ=>ob(catsub),
    :Y=>ob(frag1sym),
    :Z=>ob(frag2sym)))
  bundle_legs(out, [[1,2,3], [4,5]])
end

function multisplit_generators(enzyme::Symbol, molecule::Symbol)
  gens = Dict{Symbol, Any}()
  for ii in 1:(length(String(molecule))-1)
    gens[Symbol(:cat, enzyme, ii, sep_sym, :sub, molecule)] = enz_sub_split(enzyme, molecule, ii)  
    frag1 = Symbol(String(molecule)[1:ii])
    frag2 = Symbol(String(molecule)[ii+1:end])
    if ii != 1
      merge!(gens,multisplit_generators(enzyme, frag1))
    end
    if ii != length(String(molecule))-1
      merge!(gens,multisplit_generators(enzyme, frag2))
    end
  end
  gens
end

function gen_fragments(substrate::Symbol)
  frags = [substrate]
  for ii in 1:(length(String(substrate))-1)
    frag1 = Symbol(String(substrate)[1:ii])
    frag2 = Symbol(String(substrate)[ii+1:end])
    if ii != 1
      append!(frags,gen_fragments(frag1))
    end
    if ii != length(String(substrate))-1
      append!(frags,gen_fragments(frag2))
    end
  end
  frags
end

function multisplit_uwd(enzymes::Array{Symbol}, substrate::Symbol) 
  rel = RelationDiagram{Symbol}(0)

  substrates = gen_fragments(substrate)
  chemicals = vcat(substrates, enzymes)

  subs = add_junctions!(rel, length(substrates), variable=substrates)
  enzs = add_junctions!(rel, length(enzymes), variable=enzymes)
  nsubs = length(subs)
  nenzs = length(enzs)

  catx = add_parts!(rel, :Box, nenzs, name=[Symbol("cat$i") for i in enzymes])
  add_parts!(rel, :Port, nenzs, junction=enzs, box=catx)

  for x in 1:nenzs
    for y in 1:nenzs
      if y != x
        catxy = add_part!(rel, :Box, name=Symbol("cat$(enzymes[x])cat$(enzymes[y])"))
        add_parts!(rel, :Port, 2, junction=[enzs[x], enzs[y]], box=catxy)
      end
    end
  end

  for x in 1:nenzs
    for y in 1:nsubs
      for n in 1:(length(String(substrates[y]))-1)
        catxy = add_part!(rel, :Box, name=Symbol("cat$(enzymes[x])$(n)sub$(substrates[y])"))
        add_parts!(rel, :Port, 2, junction=[enzs[x], subs[y]], box=catxy)
      end
    end
  end
  add_parts!(rel, :OuterPort, length(chemicals), outer_junction = vcat(subs, enzs))
  rel
end




# CURRENTLY JUST A COPY OF enzyme_uwd
function multisplit_uwd(enzymes::Array{Symbol}, substrates::Array{Symbol}) # , sites::Vector{Vector{Int}}
  rel = RelationDiagram{Symbol}(0)

  chemicals = vcat(substrates, enzymes)

  subs = add_junctions!(rel, length(substrates), variable=substrates)
  enzs = add_junctions!(rel, length(enzymes), variable=enzymes)
  nsubs = length(subs)
  nenzs = length(enzs)

  catx = add_parts!(rel, :Box, nenzs, name=[Symbol("cat$i") for i in enzymes])
  add_parts!(rel, :Port, nenzs, junction=enzs, box=catx)

  for x in 1:nenzs
    for y in 1:nenzs
      if y != x
        catxy = add_part!(rel, :Box, name=Symbol("cat$(enzymes[x])cat$(enzymes[y])"))
        add_parts!(rel, :Port, 2, junction=[enzs[x], enzs[y]], box=catxy)
      end
    end
  end

  for x in 1:nenzs
    for y in 1:nsubs
      catxy = add_part!(rel, :Box, name=Symbol("cat$(enzymes[x])sub$(substrates[y])"))
      add_parts!(rel, :Port, 2, junction=[enzs[x], subs[y]], box=catxy)
    end
  end
  add_parts!(rel, :OuterPort, length(chemicals), outer_junction = vcat(subs, enzs))
  rel
end


