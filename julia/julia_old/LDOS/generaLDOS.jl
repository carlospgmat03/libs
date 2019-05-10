using ArgParse

function parse_commandline()
    cometin = ArgParseSettings()

    @add_arg_table cometin begin
        "--base", "-b"
            help = "Todo lo que esta antes del numero de miembro."
            required = true
        "--salida", "-s"
            help = "El nombre de salida"
            required = true
        "--origen", "-o"
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
nombreExterno = parsed_args["salida"]
nombreBase = parsed_args["base"]
origen = parsed_args["origen"]
delta = parsed_args["delta"]
#nombreBase = "j_q12_isi0p01_"
#nombreExterno = "jq12_isi0p01_delta0p1_sep1p0.dat"
for i in 1:20
    nombreSalida = string(nombreBase, "$i", "_todo")
    estadosCompleto = string(nombreSalida, "_estados_al_to.dat")
    fasesCompleto = string(nombreSalida, "_fases_al_to.dat")
    estadosCampo = string(nombreSalida, "_estados_ca_to.dat")
    corrida = `julia supermates.jl -S $origen -d $delta -e $estadosCompleto -f $fasesCompleto -p $estadosCampo -s $nombreExterno`
    run(corrida)
end
