Tracking bird's flight using Kalman filter
In this example, an observer tracks the trajectory defined by a bird.
We assume that the bird is following a 1D accelerated motion.

The dynamics are described by State Transition matrix as a linear system
X(t+dt) = A(dt)*X(t) + B(dt)*u(t)
Y(t) = C*X(t)

Kalmar filter formulation reference:
https://en.wikipedia.org/wiki/Kalman_filter#Example_application.2C_technical
