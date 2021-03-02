import os
import matplotlib.pyplot as plt
import numpy


def gera_figura(dados, ite, n):

    x = dados[n][ite]["ite"]

    fig, ax = plt.subplots()

    p = dados[n][ite]["L2P"]
    line1, = ax.plot(x, p, '-', linewidth=1,
                     label='Pressao')

    u = dados[n][ite]["L2U"]
    line2, = ax.plot(x, u, '--', linewidth=1,
                     label='U')

    u = dados[n][ite]["L2V"]
    line2, = ax.plot(x, u, ':', linewidth=1,
                     label='V')

    ax.legend(loc='upper right')

    ax.set_yscale("log")

    plt.savefig("./figuras/Residuo_Nx_Ny_%i_ite_sol_%i.png" %
                (n, ite), dpi=300, format="png")


def gera_figura_tempo(dados, n):

    x = dados[n].keys()
    x.sort()

    tempo = []
    for i in x:
        tempo.append(dados[n][i]["tempo"])

    if len(tempo) > 0:
        tempo = numpy.array(tempo)
        tempo[:] = tempo[:] / tempo[0]

        fig, ax = plt.subplots()

        line1, = ax.plot(x, tempo, '-', linewidth=1,
                         label='TCPU')

        ax.legend(loc='upper right')

        # ax.set_yscale("log")

        plt.savefig("./figuras/Processamento_Nx_Ny_%i.png" % (n),
                    dpi=300, format="png")


def gera_figura_tempo_unificado(dados):

    ns = dados.keys()
    ns.sort()

    tempos = {}
    iteracoes = {}

    fig, ax = plt.subplots()

    for n in ns:
        tempos[n] = []
        iteracoes[n] = []
        for i in range(50):
            try:
                tempos[n].append(dados[n][i]["tempo"])
                iteracoes[n].append(i)
            except (ValueError, KeyError):
                pass

        if len(tempos[n]) > 0:
            tempos[n] = numpy.array(tempos[n])
            tempos[n][:] = tempos[n][:] / tempos[n][0]

        ax.plot(iteracoes[n], tempos[n], '-', linewidth=1,
                label='%ix%i' % (n, n))

        ax.legend(loc='upper right')

    # ax.set_yscale("log")
    plt.savefig("./figuras/Processamento_unificado.png",
                dpi=300, format="png")


def gera_figura_iteracoes_unificado(dados):

    ns = dados.keys()
    ns.sort()

    ite_exts = {}
    iteracoes = {}

    fig, ax = plt.subplots()

    for n in ns:
        ite_exts[n] = []
        iteracoes[n] = []
        for i in range(50):
            try:
                ite_exts[n].append(float(dados[n][i]["ite"][-1]))
                iteracoes[n].append(i)
            except (ValueError, KeyError):
                pass

        #if len(ite_exts[n]) > 0:
        #    ite_exts[n] = numpy.array(ite_exts[n])
        #    ite_exts[n][:] = ite_exts[n][:] / ite_exts[n][0]

        ax.plot(iteracoes[n], ite_exts[n], '-', linewidth=1,
                label='%ix%i' % (n, n))

        ax.legend(loc='upper right')

    # ax.set_yscale("log")
    plt.savefig("./figuras/Iteracoes_unificado.png",
                dpi=300, format="png")


# ---------------------------------------
iteracoes_externas = 1000
ite_solver = list(range(11))
dados = {}

Nxs = [2**n for n in range(3, 9)]
threads = 1

for n in Nxs:
    dados[n] = {}

    for ite_s in ite_solver:

        diretorio = "./resultados/Nx_%i_Ny_%i/threads_%i/" % (n, n, threads)

        nome = (diretorio +
                "convergencia_UVP_ite_ext_"
                "%s_ite_sol_%s_Nx_Ny_%i_Threads_%i.txt" % (iteracoes_externas,
                                                           ite_s,
                                                           n, threads))

        if os.path.exists(nome):

            dados[n][ite_s] = {"ite": [],
                               "L2U": [],
                               "L2V": [],
                               "L2P": []}

            arq = open(nome)
            lines = arq.readlines()
            arq.close()

            print (nome)

            for l in lines[1:-1]:
                dados[n][ite_s]["ite"].append(int(l.split()[0]))
                dados[n][ite_s]["L2U"].append(float(l.split()[1]))
                dados[n][ite_s]["L2V"].append(float(l.split()[2]))
                dados[n][ite_s]["L2P"].append(float(l.split()[3]))
                dados[n][ite_s]["tempo"] = float(lines[-1].split()[-1])

        # gera_figura(dados, ite_s, n)
    # gera_figura_tempo(dados, n)

gera_figura_tempo_unificado(dados)
gera_figura_iteracoes_unificado(dados)
