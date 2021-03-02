import os
import subprocess
import sys

fontes = ["funcoes_abstratas.f90",
          "funcoes_alias.f90",
          'matriz_A.f90',
          'gs.f90',
          'extrapolacoes_de_u_e_v.f90',
          'fontes_subrotinas.f90',
          'class_array_subrotinas.f90',
          'residuo.f90',
          'variaveis_solvers_p.f90',
          'variaveis_solvers_u.f90',
          'variaveis_solvers_v.f90',
          'mg_gs_u.f90',
          'mg_gs_v.f90',
          'subrotinas_mg_gs.f90',
          "class_u.f90",
          "class_v.f90",
          "class_p.f90",
          "variaveis_gerais.f90",
          "navier_stokes_inout.f90",
          "experimentos_numericos.f90",
          "inicializacoes.f90",
          "NavierStokes.f90"]

executavel = "main"
if "linux" in sys.platform:
    compilador = "/usr/bin/gfortran"
    os.system("clear")
else:
    compilador = "gfortran"
    os.system("cls")

compilacao = [compilador]
for f in fontes:
    compilacao.append(f)

compilacao.append("-o")
compilacao.append(executavel)

# Openmp flag
compilacao.append("-fopenmp")

# 132 caracteres flag
compilacao.append("-ffree-line-length-none")

if "win" in sys.platform:
    executavel = executavel + ".exe"

if os.path.exists(executavel):
    os.remove(executavel)

p = subprocess.Popen(compilacao)
p.wait()

# Removendo os temporarios
for f in os.listdir("./"):
    if f.split(".")[-1] == "mod":
        os.remove(f)

# Executando e enviando a saida para o log
if os.path.exists(executavel):

    flag = ""
    if "linux" in sys.platform:
        flag = "./"

    os.system(flag + "%s" % (executavel))

else:
    print("\nExecutavel nao compilado\n")
