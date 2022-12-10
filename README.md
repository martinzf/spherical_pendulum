# spherical_pendulum
### About
MATLAB R2022a project for plotting an animated spherical pendulum given 
certain initial conditions.

### How to use
* Clone the repository.
* Open the folder from MATLAB's file explorer.
* Double click `Spherical_pendulum.prj` or open `main.m` and click `Run`.
* Answer the prompts in the `Command Window`, hitting enter after 
inputting each value.
* Wait while the animation loads, MATLAB should open a plot window.

### Preview
![Alt Text](preview.gif)

### Theory
Using generalised coordinates $\theta$ (polar angle) and $\varphi$ 
(azimuthal angle), one can compute the Lagrangian of a spherical pendulum
with a bob of mass $m$ and rod of length $l$ to be:

$$L=ml^2[\frac{1}{2}(\dot{\theta}^2+\sin^2(\theta)\dot{\varphi}^2)
-\frac{g}{l}\cos(\theta)]$$ 

This results in the following Euler-Lagrange equations:

$$\ddot{\theta}=\sin(\theta)\cos(\theta)\dot{\varphi}^2$$

$$\sin^2(\theta)\dot{\varphi}=const.\equiv c$$

Using the energy integral, the equation for $\theta$ can be simplified, 
and we obtain the equation of motion of an effective particle in a 1D
potential:

$$\frac{1}{2}\dot{\theta}^2+U(\theta)=\varepsilon$$

$$U(\theta)=\frac{c^2}{2\sin^2(\theta)}+\frac{g}{l}\cos(\theta)$$

And thus the problem is reduced to quadratures:

$$t=\pm\int\frac{d\theta}{\sqrt{\varepsilon-U(\theta)}}$$

$$\varphi=c\int\frac{dt}{\sin^2(\theta)}$$