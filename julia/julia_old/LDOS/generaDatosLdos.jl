using ArgParse

function parse_commandline()
    cometin = ArgParseSettings()

    @add_arg_table cometin begin
        "--base", "-b"
            help = "Nombre base de los archivos que tienen los estados y fases."
            required = true
        "--elementos", "-n"
            help = "Cantidad de elementos del ensamble."
            required = true
            arg_type = Int
            default = 15
        "--qbits", "-q"
            help = "NÃmero de qbits."
            required = true
            arg_type = Int
            default = 9
        "--ising", "-i"
            help = "Intensidad de la interaccion interespin."
            required = true
            arg_type = Float64
            default = 0.7
    end
    return parse_args(cometin)
end

parsed_args = parse_commandline()
nQbits = parsed_args["qbits"]
ising = parsed_args["ising"]
nombreBase = parsed_args["base"]
cElementos = parsed_args["elementos"]
for i in 1:cElementos
    nombreSalida = string(nombreBase, "$i", "_todo")
    corrida = `julia fasesyestadostodo.jl --ising $ising -q $nQbits -o $nombreSalida`
    run(corrida)
end
