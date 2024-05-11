#import Pkg
#using CondaPkg
#ENV["PYTHON"] = ""
#Pkg.build("PythonCall")
#Conda.pip_interop(true)
#Conda.pip("install", "proposal")
#add_pip("proposal"; version="")