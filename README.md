# CCBounce
This is a piece of C++ code to Calculate Bounce action and tunneling temperature in cosmological first-order phase transition scenarios.

# Requirements

It only requires a c++ compiler. 	You can adapt the following line in the makefile to your system:

"g++  -std=c++11 CCBounce.cpp -o a.out && ./a.out"

For visualization, it is required to install GNUPLOT:

MAC: brew install gnuplot

LINUX: apt-get install gnuplot

The code can be runned with the commande "make" typed in a terminal opened in the same directory.

# Usage

The software calculates the bounce action with an over-shooting/under-shooting method, using Runge-Kutta order 5 with adaptive step-size.
Then it compares it with the critical action, calculate the nucleation temperature and the phase transition completion rate \beta/H.
The error on the bounce action and nucleation temperature are calculated and below percent level.
The code embeds the minimal scale-invariant U(1) extension of the Standard Model studied in https://arxiv.org/pdf/2311.13640.pdf.
The scale-invariant SU(2)_D model is also embedded.
It is straightforward to extend the code to any other phase transition model (scale-invariant or not).
Most of the lines of code are commented to help the user to understand its structure and implement any modifications.

Without any modification of the code, you can already use it to perform four tasks (in the scale-invariant U(1) extension of the Standard Model):
1) Plot the scalar potential,
2) Calculate and plot the bouncing trajectory at a given temperature T,
3) Compute tables of bounce action for a list of values of gauge coupling constant g_X and temperature T,
4) Calculate the nucleation temperature T_n and the phase transition completion rate \beta/H.
