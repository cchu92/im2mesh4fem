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
lsh = dot(inner(sigma(du), eps(u_))*dx)
rsh = dot(f,u)*dx
```






