for.alg <- function (p.init, m.init, params_dist_emission, params_dist_sojourn, n_states, n, y) {
    ####################################################
    # Modelo con n estados latentes con distribuciones #
    # de emisión y de permanencia gamma.               #
    ####################################################
  
    # Se inicializan las probabilidades iniciales (p.init) en logaritmo
    p.init <- log(p.init)

   if (length(m.init) > 0) {
    if (all(dim(m.init) == c(n_states, n_states)) && all(diag(m.init) == 0)) {
        P <- m.init
    } else {
        stop("Error: m.init debe ser una matriz de tamaño num_states x num_states con ceros en la diagonal.")
    }
  } else {
    P <- matrix(0, n_states, n_states)
    for (i in 1:n_states) {
        if (i > 1) P[i, i-1] <- runif(1, 0, 1)  # Transición hacia atrás
        if (i < n_states) P[i, i+1] <- runif(1, 0, 1)  # Transición hacia adelante
    }
    diag(P) <- 0
  }
    P <- t(apply(P, 1, function(x) x / sum(x)))
    #print("===== Matriz de Transición =====")
    #print(P)
    # Z es una matriz de unos con la diagonal en 0, se utiliza para manejar las no-diagonales
    Z <- matrix(rep(1, n_states^2), n_states, n_states) - diag(n_states)
    
    ### Cálculo de las distribuciones de tiempos de permanencia (Gamma) ###
    log.P <- log(P)
    log.P[P == NA] <- -Inf

   # print("===== Log(Matriz de Transición) =====")
   # print(log.P)
    
   
    # Se calculan las densidades de la distribución de permanencia Gamma para cada estado
    log.h <- matrix(0, n-1, n_states)
    
    for (i in 1:n_states) {

      log.h[,i] <- dgamma(1:(n-1), shape=params_dist_sojourn[i, 1], scale=params_dist_sojourn[i, 2], log=T) - 
                   log(1 - pgamma(0, shape=params_dist_sojourn[i, 1], scale=params_dist_sojourn[i, 2]))
    }

   # print("===== Log(Densidades de distribución de permanencia) (n-1 x n_states) =====")
   # print(log.h)

    # Se inicializa la primera columna de densidades de tiempo de permanencia
    log.H1 <- log.h[1,]

   # print("===== Log(Densidades de tiempo de permanencia para T = 1) (1 x n_states) =====")
   # print(log.H1)

    # Matriz de transición multiplicada por las densidades de permanencia (para 1 paso de tiempo)
    log.PH1 <- log.P + Z %*% diag(log.H1)
    #Hacer calculos

   # print("===== Log(Matriz de transición multiplicada por las densidades de permanencia) (1 x n_states) =====")
   # print(log.PH1)

    # Cálculo de las razones entre las densidades h(t) y h(t-1) para cada estado
    Ratios <- matrix(0, n-2, n_states)
    for (i in 1:n_states) {
      Ratios[,i] <- log.h[2:(n-1),i] - log.h[1:(n-2),i]
    }


   # print("===== razones entre las densidades h(t) y h(t-1) para cada estado =====")
   # print(Ratios)

    #####

    # Cálculo de las funciones de supervivencia para la distribución Gamma (permanencia en el estado)
    log.s <- matrix(0, n, n_states)
    for (i in 1:n_states) {
      log.s[,i] <- pgamma(0:(n-1), shape=params_dist_sojourn[i, 1], scale=params_dist_sojourn[i, 2], lower.tail=F, log.p=T) - 
                   log(1 - pgamma(0, shape=params_dist_sojourn[i, 1], scale=params_dist_sojourn[i, 2]))
    }
    #####

    # Se inicializa la primera columna de supervivencia
    log.S1 <- log.s[1,]
    #####

    # Matriz de transición multiplicada por los valores de supervivencia para 1 paso de tiempo
    log.PS1 <- log.P + Z %*% diag(log.S1)
    #####

    # Cálculo de las razones entre las funciones de supervivencia y las densidades de permanencia
    S.Ratios <- matrix(0, n-1, n_states)
    for (i in 1:n_states) {
      S.Ratios[,i] <- log.s[2:n,i] - log.h[1:(n-1),i]
    }
    ####

  
    ### Cálculo de los valores del modelo condicional Gamma (emisión) ###
    log.F <- matrix(0, n, n_states)
    for (i in 1:n_states) {
      log.F[,i] <- dgamma(y, shape=params_dist_emission[i, 1], scale=params_dist_emission[i, 2], log=T)  # Densidad gamma para los valores observados
    }
    ###
  
    ### FUNCIÓN PARA CALCULAR EL LOG DE LA SUMA DE LOGS ###
    sumlog <- function(vec) {
      M <- max(vec)
      if(is.na(M)){ return(-Inf) }
      if (M == -Inf) { return(M) }
      result <- log(sum(exp(vec - M))) + M
      return(result)
    }
  
    ### INICIO DE LA RECURSIÓN ###
    alpha <- list()  # Lista para almacenar los vectores alpha
    Alpha <- list()  # Lista para almacenar las matrices alpha
    xi <- matrix(rep(0, n * n_states), n, n_states)  # Matriz para almacenar xi
    log.lik <- numeric(n)  # Vector para almacenar la verosimilitud en cada paso
  
    ### t=1 ###
    gamma <- p.init + log.H1 + log.F[1,]  # Cálculo de gamma para el primer paso
    log.lik[1] <- sumlog(gamma)  # Cálculo de la verosimilitud para t=1
    alpha[[1]] <- gamma - log.lik[1]  # Se ajusta alpha con el log-likelihood
    Alpha[[1]] <- matrix(rep(-Inf, n_states * n_states), n_states, n_states)
    diag(Alpha[[1]]) <- alpha[[1]]  

    ### t=2 ###
    gamma <- gamma + Ratios[1,] + log.F[2,]  # Cálculo de gamma para t=2 (diagonales)
    A <- diag(alpha[[1]]) %*% Z + log.PH1 + Z %*% diag(log.F[2,])  # Cálculo de A (no-diagonales)
    xi[2,] <- 0
    
    for (j in 1:n_states) {
      xi[2,j] <- sumlog(A[-j,j])
    }

    xi[2,] <- xi[2,] + log.lik[1]  # Escalado para xi
    alpha[[2]] <- apply(matrix(c(gamma, xi[2,]), n_states, 2), 1, sumlog)
    log.lik[2] <- sumlog(alpha[[2]])  # Verosimilitud para t=2
    alpha[[2]] <- alpha[[2]] - log.lik[2]  # Escalado de alpha para t=2
    Alpha[[2]] <- A + matrix(log.lik[1], n_states, n_states)
    diag(Alpha[[2]]) <- gamma
   
    ### t=3,...,n-1 ###
    for (i in 3:(n-1)) {
      gamma <- gamma + Ratios[i-1,] + log.F[i,]  # Cálculo de gamma para t>=3 (diagonales)
      A <- diag(alpha[[i-1]]) %*% Z + log.PH1 + Z %*% diag(log.F[i,])  # Cálculo de A (no-diagonales)
      
      #xi[2:(i-1),] <- t(t(xi[2:(i-1),] + Ratios[(i-2):1,]) + log.F[i,])  # Cálculo de xi (diagonales)
      
      #Revisar de nuevooo

      for (k in 2:(i-1)) {
          for(j in 1:n_states){
              xi[k, j] <- xi[k, j] + Ratios[i-k, j] + log.F[i, j]
          }
      }
      
      for (j in 1:n_states) {
        xi[i,j] <- sumlog(A[-j, j])
      }
      xi[i,] <- xi[i,] + log.lik[i-1]  # Escalado para xi

      dummy <- apply(matrix(c(gamma, xi[2,]), n_states, 2), 1, sumlog)

      for (j in 3:(i-1)) {
        dummy <- apply(matrix(c(dummy, xi[j,]), n_states, 2), 1, sumlog)
      }
      alpha[[i]] <- apply(matrix(c(dummy, xi[i,]), n_states, 2), 1, sumlog)  # Actualización de alpha
      log.lik[i] <- sumlog(alpha[[i]])  # Verosimilitud para t=i
      alpha[[i]] <- alpha[[i]] - log.lik[i]  # Escalado de alpha
      Alpha[[i]] <- A + matrix(log.lik[i-1], n_states, n_states) - matrix(log.lik[i], n_states, n_states)
      diag(Alpha[[i]]) <- dummy - log.lik[i]
    }
  
    ### t=n -- Último paso, usando las funciones de supervivencia ###
    gamma <- gamma + S.Ratios[n-1,] + log.F[n,]  # Cálculo de gamma para el último paso
    #xi[2:(n-1),] <- t(t(xi[2:(n-1),] + S.Ratios[n-2:1,]) + log.F[n,])  # Cálculo de xi para el último paso
    A <- diag(alpha[[n-1]]) %*% Z + log.PS1 + Z %*% diag(log.F[n,])  # Cálculo de A (no-diagonales)

    for (k in 2:(n-1)) {
          for(j in 1:n_states){
              xi[k, j] <- xi[k, j] + Ratios[n-k, j] + log.F[n, j]
          }
    }



    for (j in 1:n_states) {
      xi[n,j] <- sumlog(A[-j, j])
    }

    xi[n,] <- xi[n,] + log.lik[n-1]  # Escalado para xi en el último paso

    dummy <- apply(matrix(c(gamma, xi[2,]), n_states, 2), 1, sumlog)
    for (j in 3:n) {
        dummy <- apply(matrix(c(dummy, xi[j,]), n_states, 2), 1, sumlog)
    }
    alpha[[n]] <- apply(matrix(c(dummy, xi[n,]), n_states, 2), 1, sumlog)  # Actualización de alpha
    log.lik[n] <- sumlog(alpha[[n]])  # Verosimilitud total
    Alpha[[n]] <- A + matrix(log.lik[n-1], n_states, n_states) - matrix(log.lik[n], n_states, n_states)
    diag(Alpha[[n]]) <- dummy - log.lik[n]


    return(log.lik[n])
  
}
