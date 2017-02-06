using ArgParse

function parse_commandline()
    cometin = ArgParseSettings()

    @add_arg_table cometin begin
        "--estados", "-e"
            help = "Archivo que contiene los estados del sistema
            total."
            required = true
        "--fases", "-f"
            help = "Archivo que contiene las fases del sistema
            total."
            required = true
        "--patada", "-p"
            help = "Archivo que contiene los estados del sistema
            patada."
            required = true
        "--salida", "-s"
            help = "Archivo de salida"
            required = true
        "--sep", "-S"
            help = "valor de Ek"
            required = true
            arg_type = Float64
            default = 0.0
        "--delta", "-d"
            help = "Valor de delta"
            required = true
            arg_type = Float64
            default = 0.1
    end

    return parse_args(cometin)
end

parsed_args = parse_commandline()
arcEstados = parsed_args["estados"]
arcFases = parsed_args["fases"]
arcPatada = parsed_args["patada"]
arcSalida = parsed_args["salida"]
separacion = parsed_args["sep"]
delta = parsed_args["delta"]
vectoresTotal = map(x->eval(parse(x)),readcsv(arcEstados))
fasesTotal = map(x->eval(parse(x)),readcsv(arcFases))
#readdlm("./j_q14_k1_isi0p01_cal_fases_al.dat")
#diagtrucha = Diagonal(vec(exp(-1im*fasesTotal)))
diagtrucha = Diagonal(vec(fasesTotal))
mattrucha = vectoresTotal*diagtrucha*inv(vectoresTotal)
vectoresPatada = map(x->eval(parse(x)),readcsv(arcPatada))
viste = Float64[]
anglefasestotal = angle(vec(fasesTotal))
iteraciones = length(anglefasestotal)
for i in 1:iteraciones
    push!(viste, angle(dot(vectoresPatada[:,i], mattrucha*vectoresPatada[:,i])))
end
meta = (viste - mean(anglefasestotal))/sqrt(var(anglefasestotal))
metatron = (anglefasestotal - mean(anglefasestotal))/sqrt(var(anglefasestotal))
indices = Int[]
for i in 1:iteraciones
    #if abs(meta[i]) <= .025 para I <= .1 y .1 para las otras
    if abs(meta[i]-separacion) <= delta
        push!(indices, i)
    end
end
estados = Array{Complex{Float64},1}[]
for elem in indices
    push!(estados, vectoresPatada[:,elem])
end

valores = Tuple{Float64,Float64}[]
for i in 1:iteraciones
    total = 0
    for elem in estados
        total += abs(dot(vectoresTotal[:,i], elem))^2
    end
    push!(valores, (metatron[i],total))
end

fff = open(arcSalida, "a")
writedlm(fff, valores, ",")
