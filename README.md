## RUST FLUID SIMULATION / DYNAMICS

Based on Navier-Stokes equations on velocity and density of viscous fluids.
 
Heavily influenced / based-upon Jos Stam paper of the topic. Can be found
here: [www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf](www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf).
 
I'm new to Rust (this is my first time actually programming with it) so
there are bound to be tons of stylistic/conventional sins. If you've noticed
any and have the time of day, feel free to let me know about them. :)

Also I fully recognize that a lot of the calculations could be made more 
efficient by working around some of the Rust reference/memory safety 
mechanisms and avoiding some of the object-oriented / generics I used, but
I was really just using this as an opportunity to familiarize myself with
as much of Rust as I could.

In particular, I would love a way to macro out the for loops in the generic
vector structures I made, so that I am not running a 1 to 2 iteration loop 
on every... single... vector operation.

