
##' Fisher information matrix dd term with all the arguments
I_pulse <- function(pars, NB = TRUE) {
  eval(expression({
    I <- t^2 * exp(-2*d*t) * mu/(1 - exp(-d * t))
    if (NB) {
      I <- I * 1/(1 + mu/k * (1 - exp(-d * t)))
    }
    I
  }),
  envir = pars)
}

I_chase <- function(pars, NB = TRUE) {
  eval(expression({
    I <- t^2 * mu * exp(-d*t)
    if (NB) {
      I <- I /(1 + mu/k *  exp(-d*t))
    }
    I
  }), pars)
}

I_slam <- function(pars, NB = TRUE) {
  I_pulse(pars, NB = NB) +
    I_chase(pars, NB = NB)
}

I_slam_inv <- function(pars) {
  eval(expression({
  (exp(-2*d*t)*(2*mu-4*mu*exp(d*t)+
                (2*mu-k)*exp(2*d*t)+k*exp(3*d*t)))/
    (k*mu*t^2)
  }), pars)
}
