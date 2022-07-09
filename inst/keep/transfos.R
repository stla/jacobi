z <- 0.5 + 0.6i
tau <- 1 + 1i
alpha <- sqrt(-1i * tau) * exp(pi / tau * 1i * z * z)

jacobi:::jtheta1_cpp(z, tau)
1i * jacobi:::jtheta1_cpp(z / tau, -1 / tau) / alpha 

jacobi:::jtheta2_cpp(z, tau)
jacobi:::jtheta4_cpp(z / tau, -1 / tau) / alpha

jacobi:::jtheta3_cpp(z, tau)
jacobi:::jtheta3_cpp(z / tau, -1 / tau) / alpha

jacobi:::jtheta4_cpp(z, tau)
jacobi:::jtheta2_cpp(z / tau, -1 / tau) / alpha

z <- 1 - 1i
tau <- 1e-10i
alpha <- sqrt(-1i * tau) * exp(pi / tau * 1i * z * z / pi / pi)
jacobi:::altjtheta1(z, tau)
1i * jacobi:::altjtheta1(z / tau, -1 / tau) / alpha 
jacobi::jtheta1(z, tau)
1i * jacobi::jtheta1(z / tau, -1 / tau) / alpha 


