### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 86ffd357-1510-4d05-8a38-b59b42b79b39
begin
using Pkg
Pkg.activate(".")
end

# ╔═╡ 3779b846-e5ec-4239-a1d4-af2f8c2f10eb
begin
  include("EnzymeReactions.jl")
  include("AffinityViz.jl")
  import .EnzymeReactions: ob, ode,
               inactivate, bindunbind, degrade,
               enzX, enzXY, enzXsubY,
               enz, enz_enz, enz_sub,
               enzyme_uwd
  using PlutoUI

  using AlgebraicPetri
  using Catlab.WiringDiagrams
  using Catlab.CategoricalAlgebra
  using Catlab.Graphics
  using LabelledArrays

  using DifferentialEquations
  using Plots
  using Interact
  using JSON2
  display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>".75"));
  nothing
end

# ╔═╡ 93df89f0-8429-4fcc-bd01-6982417f5134
begin
  # Initial Concentrations
  K = :K=>33000;
  S = :S=>33000;
  L = :L=>33000;
  Kinact = :K_inact=>0;
  Sinact = :S_inact=>0;
  Linact = :L_inact=>0;
  E = :E=>700000;
  G = :G=>1300000;

  # Parameter Rates (units of pM and min)
  rxns = Dict(
    :K => [inactivate(K, 7.494e-10)
           bindunbind(K, K, 7.814e-4, 3.867e-3)
           degrade(K, K, 2.265e-1)
           bindunbind(K, Kinact, 7.814e-4, 3.867e-3)
           degrade(K, Kinact, 2.265e-1)],
    :S => [inactivate(S, 1.906e-2)
           bindunbind(S, S, 3.534e-7, 1.688e2)
           degrade(S, S, 1.433e1)
           bindunbind(S, Sinact, 3.534e-7, 1.688e2)
           degrade(S, Sinact, 1.433e1)],
    :L => [inactivate(L, 7.810e-3)
           bindunbind(L, L, 1.000e-11, 7.440e3)
           degrade(L, L, 2.670e2)
           bindunbind(L, Linact, 1.000e-11, 7.440e3)
           degrade(L, Linact, 2.670e2)],
    :KE => [bindunbind(K, E, 9.668e-6,1.000e-2)
            degrade(K, E, 1.728e0)],
    :KG => [bindunbind(K, G, 2.764e-6, 8.780e-1)
            degrade(K, G, 1.502e0)],
    :SE => [bindunbind(S, E, 4.197e-7, 1.06e-3)
            degrade(S, E, 1.384e4)],
    :SG => [bindunbind(S, G, 5.152e-8, 3.894e-3)
            degrade(S, G, 8.755e-1)],
    :LE => [bindunbind(L, E, 1.977e-8, 1.000e-2)
            degrade(L, E, 1.066e2)],
    :LG => [bindunbind(L, G, 3.394e-8, 2.365e1)
            degrade(L, G, 4.352e0)],
    :KS => [bindunbind(K, S, 8.822e-4, 4.114e5)
            degrade(K, S, 9.000e-10)
            bindunbind(K, Sinact, 8.822e-4, 4.114e5)
            degrade(K, Sinact, 9.000e-10)],
    :KL => [bindunbind(K, L, 1.756e-4, 3.729e4)
            degrade(K, L, 6.505e6)
            bindunbind(K, Linact, 1.756e-4, 3.729e4)
            degrade(K, Linact, 6.505e6)],
    :SK => [bindunbind(S, K, 3.679e-4, 1.562e3)
            degrade(S, K, 4.410e2)
            bindunbind(S, Kinact, 3.679e-4, 1.562e3)
            degrade(S, Kinact, 4.410e2)],
    :SL => [bindunbind(S, L, 1.000e-3, 5.000e2)
            degrade(S, L, 1.000e-7)
            bindunbind(S, Linact, 1.000e-3, 5.000e2)
            degrade(S, Linact, 1.000e-7)],
    :LK => [bindunbind(L, K, 1.000e-3, 4.118e3)
            degrade(L, K, 3.234e1)
            bindunbind(L, Kinact, 1.000e-3, 4.118e3)
            degrade(L, Kinact, 3.234e1)],
    :LS => [bindunbind(L, S, 1.056e-12, 5.000e2)
            degrade(L, S, 5.000e-1)
            bindunbind(L, Sinact, 1.056e-12, 5.000e2)
            degrade(L, Sinact, 5.000e-1)]
  );

  # define labels to reaction network mappings
  functor(x) = oapply(x, Dict(
    :catK=>enz(rxns, K),
    :catS=>enz(rxns, S),
    :catL=>enz(rxns, L),
    :catKcatS=>enz_enz(rxns, K,S),
    :catKcatL=>enz_enz(rxns, K,L),
    :catScatK=>enz_enz(rxns, S,K),
    :catScatL=>enz_enz(rxns, S,L),
    :catLcatK=>enz_enz(rxns, L,K),
    :catLcatS=>enz_enz(rxns, L,S),
    :catKsubE=>enz_sub(rxns, K,E),
    :catSsubE=>enz_sub(rxns, S,E),
    :catLsubE=>enz_sub(rxns, L,E),
    :catKsubG=>enz_sub(rxns, K,G),
    :catSsubG=>enz_sub(rxns, S,G),
    :catLsubG=>enz_sub(rxns, L,G)));
  nothing
end

# ╔═╡ e6589d31-dce7-42c3-b494-db03fe561ae9
  uwd = enzyme_uwd([:K], [:G]);

# ╔═╡ 7dbe9349-8b9e-4ac2-b4bf-b59f58a10ebc
begin
  display_uwd(uwd)
end

# ╔═╡ 8509bb5c-84bc-43fd-a4c6-5801b210bf58
begin 
	function prepData(a,b,c,e)
		d = Dict("rates" =>[
			Dict([("slider", 1), ("max", 1), ("min",0), ("sname","rate1")]),
			],
	"conc" =>[Dict([("slider", 1), ("max", 1), ("min",0), ("sname","c1")]),
			],
		"tab" => "null",
		)
		println(d)
		
		for i = 1:length(a)
			if i ==1
				d["rates"][1] = Dict([("slider", b[1]), ("max", 1), ("min",0), ("sname",string(a[i]))]) 
				println(d)
			else
				push!(d["rates"],Dict([("slider", b[1]), ("max", 1), ("min",0), ("sname",string(a[i]))]))
			end
		end
			
		for i = 1:length(c)
			if i ==1
				d["conc"][1] = Dict([("slider", e[i]), ("max", 1), ("min",0), ("sname",string(c[i]))]) 
			else
				push!(d["conc"],Dict([("slider", e[i]), ("max", 1), ("min",0), ("sname",string(c[i]))]))
			end
		end
		
		return JSON2.write(d)
	end
end

# ╔═╡ 28ff5981-52b8-40d5-9439-98acf84018e7
begin
	function sliderInterface(dat)
		println(dat)
		return HTML("""
<html>
<head>
  <link href="https://fonts.googleapis.com/css?family=Roboto:100,300,400,500,700,900" rel="stylesheet">
  <link href="https://cdn.jsdelivr.net/npm/@mdi/font@4.x/css/materialdesignicons.min.css" rel="stylesheet">
  <link href="https://cdn.jsdelivr.net/npm/vuetify@2.x/dist/vuetify.min.css" rel="stylesheet">
  <meta name="viewport" content="width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no, minimal-ui">
</head>
<body>
  <div id="app">
    <v-app>
      <v-main>
        <template>
          <v-card>
            <v-tabs
              v-model="tab"
              background-color="transparent"

              grow
            >
              <v-tab>
                Rates
              </v-tab>
              <v-tab>
                Concentrations
              </v-tab>

            </v-tabs>

            <v-tabs-items v-model="tab">
              <v-tab-item
              >
                <v-card  flat>
                  <v-card-text>Adjust the reaction rates</v-card-text>

                  <v-slider
                    v-for = "s in rates"
                    v-model="s.slider"
                    :max="s.max"
                    :min="s.min"
					step =".00001"
                    :label = "s.sname"
                  >
                    <template v-slot:append>
                      <v-text-field
                        v-model="s.slider"
                        class="mt-0 pt-0"
                        hide-details
                        single-line
                        type="number"
                        style="width: 60px"
                      ></v-text-field>
                    </template>
                </v-slider>

    </v-slider>
  </v-card>
</v-tab-item>

<v-tab-item
>
  <v-card  flat>
    <v-card-text>Adjust the Concentrations</v-card-text>

    <v-slider
      v-for = "s in conc"
      v-model="s.slider"
      :max="s.max"
      :min="s.min"
      step =".00001"
      :label = "s.sname"
    >
      <template v-slot:append>
        <v-text-field
          v-model="s.slider"
          class="mt-0 pt-0"
          hide-details
          single-line
          type="number"
          style="width: 60px"
        ></v-text-field>
      </template>
  </v-slider>

</v-slider>
</v-card>
</v-tab-item>


</v-tabs-items>
</v-card>
</template>
<v-btn>Update </v-btn>
      </v-main>
    </v-app>
  </div>

  <script src="https://cdn.jsdelivr.net/npm/vue@2.x/dist/vue.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/vuetify@2.x/dist/vuetify.js"></script>
  <script>
    new Vue({
      el: '#app',
      vuetify: new Vuetify(),
      data () {
      return $dat
			},
    })
  </script>
</body>
</html>

			
			""")
	end
end

# ╔═╡ ee9159e2-67b6-4503-a338-10337144b821
# sliderInterface(dataPrepped)

# ╔═╡ cf9e03db-42b7-41f6-80ce-4b12ddb93211
begin
  model = uwd |> functor |> apex;
  r = join(["'$k': $(rates(model)[k])" for k in keys(rates(model))], ", ");
	
  r2 = join(["'$k': $(concentrations(model)[k])" for k in keys(concentrations(model))], ", ");
  #== Basic Slider Template
<li><input
  class='slider'
  type='range'
  min='-10'
  step='0.01'
  max='10'
  value='0'
  oninput='this.nextElementSibling.value=(10**this.value).toExponential(2)'
  ></input>
<output>1</output></li>
==#

form_vals = HTML("""
<form>
		
  <div id ="myDIV" style='height:500px;overflow:scroll'>
		
  <h4>Adjust Rates</h4>
  <table>
  </table>
  


  </div>
 <div id ="concDiv" style='height:500px;overflow:scroll'>
		<h4>Adjust Concentrations</h4>
		<table id = "table2">
		</table>
	</div>

 <div id ="bindingAff" style='height:500px;overflow:scroll'>
		<h4>Binding Affinities</h4>
		<table id = "bindingAffTable">
		</table>
	</div>
		
		
</form>

<style>
		
#myDIV {
  width: 100%;
  padding: 50px 0;
  text-align: center;
  background-color: #B3A369;
  margin-top: 20px;
}
#concDiv {
  width: 100%;
  padding: 50px 0;
  text-align: center;
  background-color:  #B3A369;
  margin-top: 20px;
  display: none
}

#bindingAff {
  width: 100%;
  padding: 50px 0;
  text-align: center;
  background-color:  #B3A369;
  margin-top: 20px;
  display: none
}
		
.button {
  background-color: #003057;
  border: none;
  color: white;
  padding: 15px 32px;
  text-align: center;
  text-decoration: none;
  display: inline-block;
  font-size: 16px;
  margin: 4px 2px;
  cursor: pointer;
}
.button:focus {     
    background-color:#BFB37C;    
}

.slider {
  -webkit-appearance: none;
  width: 100%;
  height: 7px;
  border-radius: 5px;  
  background: #d3d3d3;
  outline: none;
  opacity: 0.7;
  -webkit-transition: .2s;
  transition: opacity .2s;
}

.slider::-webkit-slider-thumb {
  -webkit-appearance: none;
  appearance: none;
  width: 20px;
  height: 20px;
  border-radius: 50%; 
  background: #003057;
  cursor: pointer;
}

.slider::-moz-range-thumb {
  width: 20px;
  height: 20px;
  border-radius: 50%;
  background: #003057;
  cursor: pointer;
}
		
</style>
		

		
<script>
//`currentScript` is the current script tag - we use it to select elements//

function hideShowRate() {
	var x = document.getElementById("myDIV");
	var y = document.getElementById("concDiv");
	var z = document.getElementById("bindingAff");
	z.style.display = "none";
	x.style.display = "block";
	y.style.display = "none"
		
}

function hideShowConc() {
	var x = document.getElementById("concDiv");
	var y = document.getElementById("myDIV");
	var z = document.getElementById("bindingAff");
	x.style.display = "block";
	y.style.display = "none"
	z.style.display = "none"

}

function hideShowBind() {
	var x = document.getElementById("bindingAff");
	var y = document.getElementById("myDIV");
	var z = document.getElementById("concDiv");
	x.style.display = "block";
	y.style.display = "none"
	z.style.display = "none"
	generateBindingTable()

}
		
function generateBindingTable(){
		
		let keyNames = Object.keys(rates)
		let objVals = Object.values(rates)
		let bindArr = []
		let unbindArr = []
		
		for (let i in keyNames){
		   
		   let nm = keyNames[i]
			console.log(nm.substring(0,4))
			if(nm.substring(0,4) == 'bind'){
				bindArr.push(nm)
		}else if(nm.substring(0,6) == 'unbind'){
				unbindArr.push(nm)
			}
		    
			
			}

		let kdArr = []; //array of binding affinities
		let kdNames = []; // array of binding affinity names
			//grab values from sliders
	    let x = form.getElementsByClassName('sliderval');
		let xarr = Array.from(x, (v,_)=>{return v.value})
		
		//Loop over bind and unbind arrays to do calculations
		for(let i in keyNames){
			
			if(keyNames[i].substring(0,4) == 'bind'){
		console.log(keyNames[i].substring(0,4))
			let bind = keyNames[i]
			for(let j in keyNames){
		    if(keyNames[j].substring(0,6) == 'unbind'){
		console.log(keyNames[j].substring(0,4))
				let unbind = keyNames[j]
					if(bind.substring(5,bind.length) == unbind.substring(7,unbind.length)){
						kdArr.push(xarr[i]/xarr[j])
						kdNames.push('Kd_'+bind.substring(5,bind.length))
					}}
					}
				}
			}

		console.log(kdArr)
		console.log(kdNames)
		const list3 = form.querySelector('#bindingAffTable')
		while(list3.rows.length > 0) {
				  list3.deleteRow(0);
				}


		for ( var r in kdNames ){
		  
		  var item = document.createElement('tr');
		  var label = document.createElement('th');
		  label.innerText = kdNames[r] +': ' + kdArr[r]

		  item.appendChild(label)

		  list3.appendChild(item)
		}
		
		}
const form = currentScript.parentElement.querySelector('form')
const list = form.querySelector('table')
var rates = {$r};
var conc = {$r2};


		
		
for ( var r in rates ){
  console.log(r);
  var item = document.createElement('tr');
  var label = document.createElement('th');
  
  label.innerText = r
  var slider_box = document.createElement('td');
  var slider_val_box = document.createElement('td');
  var slider = document.createElement('input');
  var slider_val = document.createElement('input');
  slider_val.setAttribute('type', 'text');
  slider_val.setAttribute('class', 'sliderval');
  slider_val.value = rates[r].toExponential(2);
  slider.setAttribute('type', 'range');
  slider.setAttribute('class', 'slider');
  slider.setAttribute('min', '-9.99');
  slider.setAttribute('max', '10');
  slider.setAttribute('value', '0.0');
  slider.setAttribute('step', '0.01');
  slider.setAttribute('oninput', `this.parentElement.nextElementSibling.children[0].value=((10**this.value)*\${rates[r]}).toExponential(2)`);
  slider_box.appendChild(slider)
  slider_val_box.appendChild(slider_val)
  item.appendChild(label)
  item.appendChild(slider_box)
  item.appendChild(slider_val_box)
  list.appendChild(item)
}

const list2 = form.querySelector('#table2')
for ( var r in conc ){
  console.log(r);
  var item = document.createElement('tr');
  var label = document.createElement('th');
  label.innerText = r
  var slider_box = document.createElement('td');
  var slider_val_box = document.createElement('td');
  var slider = document.createElement('input');
  var slider_val = document.createElement('input');
  slider_val.setAttribute('type', 'text');
  slider_val.setAttribute('class', 'sliderval');
  slider_val.value = conc[r].toExponential(2);
  slider.setAttribute('type', 'range');
  slider.setAttribute('class', 'slider');
  slider.setAttribute('min', '-9.99');
  slider.setAttribute('max', '10');
  slider.setAttribute('value', '0.0');
  slider.setAttribute('step', '0.01');
  slider.setAttribute('oninput', `this.parentElement.nextElementSibling.children[0].value=((10**this.value)*\${conc[r]}).toExponential(2)`);
  slider_box.appendChild(slider)
  slider_val_box.appendChild(slider_val)
  item.appendChild(label)
  item.appendChild(slider_box)
  item.appendChild(slider_val_box)
  list2.appendChild(item)
}

		
		
var x = form.getElementsByClassName('sliderval');
function onsubmit(){
  // We send the value back to Julia //
  form.value = Array.from(x, (v,_)=>{return v.value})
  form.dispatchEvent(new CustomEvent('input'))
  console.log(form.value)
}
var b = document.createElement('input');
b.setAttribute('type', 'button');
b.setAttribute('class', 'button');
b.value = 'Update Plot';
b.addEventListener('click', function() {onsubmit();console.log('hello from button')})
form.appendChild(b)
onsubmit()

let sp2 = document.getElementById("myDIV")
		
var d = document.createElement('input')
d.setAttribute('type','button')
d.setAttribute('class','button')
d.value = 'Rates'
d.addEventListener('click', function() {hideShowRate();})
form.insertBefore(d,sp2)
		
var e = document.createElement('input')
e.setAttribute('type','button')
e.setAttribute('class','button')
e.value = 'Concentrations'
e.addEventListener('click', function() {hideShowConc();})
form.insertBefore(e,sp2)
		
var f = document.createElement('input')
f.setAttribute('type','button')
f.setAttribute('class','button')
f.value = 'Binding Affinities'
f.addEventListener('click', function() {hideShowBind();})
form.insertBefore(f,sp2)		

</script>
""");
nothing
end

# ╔═╡ 4245e6c8-d35f-4e04-b5ea-96602ceeee3e
begin
a = []
b = []
carr = []
earr = []
for k in keys(rates(model))
	# println(k)
	push!(a,k)
	push!(b,rates(model)[k])
	# println(rates(model)[k])
end
	
for k in keys(concentrations(model))
	# println(k)
	push!(carr,k)
	push!(earr,concentrations(model)[k])
	# println(rates(model)[k])
end

	dataPrepped = prepData(a,b,carr,earr);
end

# ╔═╡ ba87cd7e-e9c7-4a20-99be-eee794f968a1
@bind c form_vals

# ╔═╡ 066b7505-e21b-467e-86c1-cea1ff80246e
begin
  cur_rate = Dict(tnames(model)[i]=>parse(Float64, c[i]) for i in 1:length(tnames(model)))
  # cur_conc = concentrations(model)
  cur_conc_prep =  Dict(snames(model)[i]=>parse(Float64, c[i+length(tnames(model))]) for i in 1:length(snames(model)))
	
	keysa = keys(cur_conc_prep)
	keyArr = Symbol[]
	valArr = Float64[]
	for key in keys(cur_conc_prep)
		push!(keyArr,key)
		push!(valArr,cur_conc_prep[key])
		
	end

  cur_conc = @LArray valArr Tuple(keyArr)
  vf = vectorfield(model);
  nothing
end

# ╔═╡ 1ba7bbe5-7a85-454e-a9cf-deaf5f00d6ad
sol = solve(ODEProblem(vf, cur_conc, (0.0,120.0),cur_rate));

# ╔═╡ a141cd27-6ea0-4f73-80b5-72d8e5770ed4
begin
  tsteps = range(0.0,120.0, length=5000)
  labels = [:G_deg]
  plot([[sol(t)[l] for t in tsteps] for l in labels], labels=hcat(String.(labels)...), linewidth=3)
end

# ╔═╡ 9625798a-67df-49e4-91ce-c7e23ed2a177
model |> AffinityNet |> to_graphviz

# ╔═╡ Cell order:
# ╠═86ffd357-1510-4d05-8a38-b59b42b79b39
# ╠═3779b846-e5ec-4239-a1d4-af2f8c2f10eb
# ╟─93df89f0-8429-4fcc-bd01-6982417f5134
# ╠═e6589d31-dce7-42c3-b494-db03fe561ae9
# ╟─7dbe9349-8b9e-4ac2-b4bf-b59f58a10ebc
# ╟─4245e6c8-d35f-4e04-b5ea-96602ceeee3e
# ╟─8509bb5c-84bc-43fd-a4c6-5801b210bf58
# ╟─28ff5981-52b8-40d5-9439-98acf84018e7
# ╠═ee9159e2-67b6-4503-a338-10337144b821
# ╠═cf9e03db-42b7-41f6-80ce-4b12ddb93211
# ╟─ba87cd7e-e9c7-4a20-99be-eee794f968a1
# ╟─066b7505-e21b-467e-86c1-cea1ff80246e
# ╟─1ba7bbe5-7a85-454e-a9cf-deaf5f00d6ad
# ╟─a141cd27-6ea0-4f73-80b5-72d8e5770ed4
# ╠═9625798a-67df-49e4-91ce-c7e23ed2a177
