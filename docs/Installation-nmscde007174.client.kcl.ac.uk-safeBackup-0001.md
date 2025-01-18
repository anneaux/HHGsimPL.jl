To use ATIsimPL.jl and/or HHGsimPL.jl you'll need to install the programming language julia and make ir available for jupyter notebooks to use (aka installing a julia kernel).

## Using julia in jupyter notebooks
- install the julia language as described [here](https://julialang.org/downloads/):
    - for Linux: `curl -fsSL https://install.julialang.org | sh` 
    - for Windows it seems you can just download it via the Microsoft Store [here](https://apps.microsoft.com/detail/9njnww8pvkmn?amp%3Bgl=gb&hl=en-GB&gl=GB) (disclaimer: I never used Windows...)
- make sure you can use jupyter notebook installed somehow
- install a julia kernel to be used by jupyter notebooks, described e.g. [here:](https://www.kdnuggets.com/2022/11/setup-julia-jupyter-notebook.html)
    1) Open a terminal (Linux) / powershell (? Windows) and enter `julia` (you should see the logo coming up and a julia prompt. This is called the julia REPL mode - as opposed to having a file "test.jl" and running it e.g. from VS code or via the terminal using `julia test.jl`)
    2) Enter the julia "package mode" by entering closed square bracket (this key: `]` ), and add the package that provides the julia kernel for jupyter by entering `add IJulia`. Alternatively, enter `using Pkg; Pkg.add("IJulia")` directly in the julia prompt.
    3) Optional: if you run code that would benefit from being parallelised (e.g. a `for` loop where the iteration steps are independent of each other) it makes sense to install a julia kernel that uses multiple cores ("threads"). For that, first, find out how many processors your laptop has (on Linux: `lscpu`, CPU, on Windows: check system settings (?)). Then, in the julia REPL, do: 
    ```
    using IJulia
    installkernel("Julia (n threads)", env=Dict("JULIA_NUM_THREADS"=>"$n"))
    ```
    (and replace `n` with a reasonable number somewhere around half of the number of your CPUs).
- open an empty jupyter notebook, you should now be able to select a julia kernel. Test e.g. by writing `println("Hello, world")` and running that cell.

## Running ATIsimPL / HHGsimPL
- get the code from the github repo I sent you: Select the `user` branch in the drop-down menu (where it says `master`) and get the code by either cloning the repo (if you know git) or simply download the code as a zip and unpack it.
- Because I don't know how to organise all the dependencies between the libraries properly, and ultimately, because my code isn't a proper package itself (yet), you'll have to add all the libraries (packages) that are used manually.
Do so by opening the julia REPL, enter package mode and then add all the following packages:
```
add StaticArrays, NLsolve, Contour, ...
```
- open the test notebook and see if it works.

