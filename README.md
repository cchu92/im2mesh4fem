inspired by the [im2mesh](https://fr.mathworks.com/matlabcentral/fileexchange/71772-im2mesh-2d-image-to-triangular-meshes) matlab code.

# Images to mesh 


# FEM simulation


# FEnics on win

Fenics is primary developed for linux and macos, for win, insall use docer

1. install docker 
2. install fenics iamges 'docker pull quay.io/fenicsproject/stable:current'
3. Mounting Local Directories 'path' and Jupyter notebook:

  docker run -ti --rm -p 8888:8888 -v llocal_directories_to_mount:/home/fenics/shared quay.io/fenicsproject/stable:current 'jupyter-notebook --ip=0.0.0.0 --port=8888 --no-browser --allow-root' 
4. lauch jupyter notebook:  http://127.0.0.1:8888/?token=your_token


# Arc length solution for D-F
[arc legnth solution](https://github.com/pprachas/fenics_arclength)



# FEniCs math 

## linear problem


$\text{Find } \boldsymbol{u}\in V \text{ s.t. } \int_{\Omega}
\boldsymbol{\sigma}(\boldsymbol{u}):\boldsymbol{\varepsilon}(\boldsymbol{v}) d\Omega
= \int_{\Omega} \boldsymbol{f}\cdot\boldsymbol{v}  d\Omega \quad \forall\boldsymbol{v} \in V
$
```
V = VectorFunctionSpace(mesh, "Lagrange", 1)
u = TrialFunction(V)
v = TestFunction(V)
lsh = dot(inner(sigma(u), eps(v))*dx)
rsh = dot(f,u)*dx
```



## weak coupled thmermomechancial problem

\begin{equation}
\boldsymbol{\sigma} = \mathbb{C}:(\boldsymbol{\varepsilon}-\alpha(T-T_0)\boldsymbol{1}) = \lambda\text{tr}(\boldsymbol{\varepsilon})\boldsymbol{1}+2\mu\boldsymbol{\varepsilon} -\alpha(3\lambda+2\mu)(T-T_0)\boldsymbol{1}
\end{equation}


```
V = VectorFunctionSpace(mesh, 'CG', 2)
u = TrialFunction(V)
v = TestFunction(V)
Wint = inner(sigma(u, Delta_T), eps(v))*dx
aM = lhs(Wint)
LM = rhs(Wint) + inner(f, v)*dx
bcu = DirichletBC(V, Constant((0., 0.)), lateral_sides)
u = Function(V, name="Displacement")
solve(aM == LM, u, bcu)
```

