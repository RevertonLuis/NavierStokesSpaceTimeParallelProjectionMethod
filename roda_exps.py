import os
import shutil

arq = open("configuracoes_gerais.txt.bkp")
lines = arq.readlines()
arq.close()

dir_main = "./resultados/"

ite_sol = list(range(1, 12))
Nxs = [2**n for n in range(3, 9)]

Nxs = [8, 16, 32, 64, 128, 256, 512, 1024, 2048]

tn = 9
iteracoes_externas = 1000

threads = [1]

for thread in threads:
    for nx in Nxs:
        for ite in ite_sol:
            arq = open("configuracoes_gerais.txt", 'w')

            for l in lines:
                l_new = l.split()
                if "Iteracoes_Solver" in l_new:
                    l_new[-1] = str(ite)

                if "Iteracoes_Externas" in l_new:
                    l_new[-1] = str(iteracoes_externas)

                if len(l_new) == 2 and "Nx" in l_new:
                    l_new[-1] = str(nx)

                if len(l_new) == 2 and "Ny" in l_new:
                    l_new[-1] = str(nx)

                if len(l_new) == 2 and "Threads" in l_new:
                    l_new[-1] = str(thread)

                if len(l_new) == 2 and "Tn" in l_new:
                    l_new[-1] = str(tn)

                arq.write(" ".join(l_new) + "\n")
            arq.close()
            # os.system("./main")
            os.system("main.exe")

            dir_nx = dir_main + "Nx_%i_Ny_%i/" % (nx, nx)
            dir_threads = dir_nx + "threads_%i/" % thread
            if not os.path.exists(dir_threads):
                os.makedirs(dir_threads)

            for a in os.listdir(dir_main):
                if a.split(".")[-1] == "txt":
                    shutil.move(dir_main + a, dir_threads + a)
