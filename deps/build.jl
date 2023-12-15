using Conda

ENV["PYTHON"] = ""
Pkg.build("PyCall")
Conda.pip_interop(true)
Conda.pip("install", "proposal")