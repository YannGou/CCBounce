## CCBounce: A C++ code to Calculate Bounce action and tunneling temperature in cosmological first-order phase transition scenarios.
[**Requirements**](#Requirements)
| [**Usage**](#Usage)
| [**Cite**](#Cite)

## Requirements

It only requires a c++ compiler. 	You can adapt the following line in the makefile to your system:

```console
$ g++  -std=c++11 CCBounce.cpp -o a.out && ./a.out
```

For visualization, it is required to install GNUPLOT:

MAC: 
```console
$ brew install gnuplot
```

LINUX:
```console
$ apt-get install gnuplot
```

The code can be simply runned with the commande "make" typed in a terminal opened in the same directory.

## Usage

The software calculates the bounce action with an over-shooting/under-shooting method, using Runge-Kutta order 5 with adaptive step-size.
Then it compares it with the critical action, calculate the nucleation temperature and the phase transition completion rate \beta/H.
The error on the bounce action and nucleation temperature are calculated and below percent level.
The code embeds the minimal scale-invariant U(1) extension of the Standard Model studied in https://arxiv.org/pdf/2311.13640.pdf.
The scale-invariant SU(2)_D model is also embedded.
It is straightforward to extend the code to any other phase transition model (scale-invariant or not).
Most of the lines of code are commented to help the user to understand its structure and implement any modifications.

Without any modification of the code, you can already use it to perform four tasks (in the scale-invariant U(1) extension of the Standard Model):
1) Plot the scalar potential,


<figure>
  <img src="https://github.com/YannGou/CCBounce/blob/d9736666eb098b4879c6006cb54d5e37a01658be/output/bounce_action.png" width="500" align="center">
  <figcaption align="center">
  Coleman-Weinberg potential
  </figcaption>
</figure>
<br/><br/>

[https://www.dropbox.com/s/e4l6kx6nhxwxdeh/bounce_action.eps?dl=0](https://www.dropbox.com/scl/fi/h6lybs59lrlbd08430h5l/potential.eps?rlkey=lo0u0v3zmp5r5y3a79ii8vsuh&dl=0)


2) Calculate and plot the tunneling trajectory at a given temperature T,

[https://www.dropbox.com/s/8ttyobok063cd3o/bounce_traj.eps?dl=0](https://www.dropbox.com/s/8ttyobok063cd3o/bounce_traj.eps?dl=0)

   
3) Compute tables of bounce action values for a list of gauge coupling constant g_X and temperature T values,
4) Calculate the nucleation temperature T_n and the phase transition completion rate \beta/H.

[https://www.dropbox.com/scl/fi/vm40qyy4acy7r617zhptp/bounce_traj.eps?rlkey=90qt0fqhf3lomy9fzpxdrlia8&dl=0](https://www.dropbox.com/scl/fi/7m346e8642vtdl10efhqv/bounce_action.eps?rlkey=vrdmx7jsh5tzgyb1sf1ylp6bt&dl=0)

**This is a research project. Expect bugs, report bugs, fix more bugs than you
create.**

## Cite
Please cite with bibtex:

@article{Gouttenoire:2023pxh,
    author = "Gouttenoire, Yann",
    title = "{Primordial Black Holes from Conformal Higgs}",
    eprint = "2311.13640",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    month = "11",
    year = "2023"
}
