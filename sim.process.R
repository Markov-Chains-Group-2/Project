sim.process <- function(L, id, num_states, noise_sd = 5, p.init, m.init, params_dist_emission, params_dist_sojourn) {
  ###############################################
  # Simula un proceso oculto de semi-Markov con #
  # transiciones entre estados según la matriz  #
  # de transición del modelo, y asigna valores  #
  # de glucosa y duraciones basados en          #
  # distribuciones específicas para cada estado.#
  ###############################################
  
  # Inicializamos las estructuras para guardar resultados
  chain <- c()    # Almacenará la secuencia de estados latentes
  time <- c()     # Almacenará los tiempos simulados
  gl <- c()       # Almacenará los valores simulados de glucosa
  
  # Definir la distribución inicial
    if (length(p.init) > 0) {
        p.init <- p.init
    } else {
        p.init <- rep(1/num_states, num_states)
    }
  
  # Matriz de transición (las transiciones permitidas son a estados adyacentes)
  if (length(m.init) > 0) {
    if (all(dim(m.init) == c(num_states, num_states)) && all(diag(m.init) == 0)) {
        # Asignar P a m.init si cumple las condiciones
        P <- m.init
    } else {
        stop("Error: m.init debe ser una matriz de tamaño num_states x num_states con ceros en la diagonal.")
    }
  } else {
    # Si m.init está vacío, crear la matriz aleatoria
    P <- matrix(0, num_states, num_states)
    for (i in 1:num_states) {
        if (i > 1) P[i, i-1] <- runif(1, 0, 1)  # Transición hacia atrás
        if (i < num_states) P[i, i+1] <- runif(1, 0, 1)  # Transición hacia adelante
    }
  }

  # Normalización de la matriz de transición
  P <- t(apply(P, 1, function(x) x / sum(x)))
  
  # Extraer los parámetros para las distribuciones de emisión desde params_dist_emission
  gamma_shape_gl <- sapply(1:num_states, function(i) params_dist_emission[i, 1])
  gamma_scale_gl <- sapply(1:num_states, function(i) params_dist_emission[i, 2])
  
  # Extraer los parámetros para las distribuciones de permanencia desde params_dist_sojourn
  gamma_shape_dur <- sapply(1:num_states, function(i) params_dist_sojourn[i, 1])
  gamma_scale_dur <- sapply(1:num_states, function(i) params_dist_sojourn[i, 2])
  
  current_state <- sample(1:num_states, 1, prob = p.init)
  
  start_time <- as.POSIXct(Sys.Date()) - sample(1:3650, 1) * 24 * 60 * 60
  start_time <- start_time + sample(0:(24*60*60 - 1), 1)
  
  while (length(chain) < L) {

    duration <- round(rgamma(1, shape = gamma_shape_dur[current_state], scale = gamma_scale_dur[current_state]))
    duration <- min(duration, L - length(chain))
    
    chain <- c(chain, rep(current_state, duration))
    
    gl_values <- rgamma(duration, shape = gamma_shape_gl[current_state], scale = gamma_scale_gl[current_state])
    gl <- c(gl, gl_values)
    
    # Generamos las marcas de tiempo en formato requerido
    if (length(time) == 0) {
      time <- seq(from = start_time, by = 5 * 60, length.out = duration)
    } else {
      last_time <- tail(time, 1)
      time <- c(time, seq(from = last_time + 5 * 60, by = 5 * 60, length.out = duration))
    }
    
    # Elegir el siguiente estado basado en las probabilidades de transición
    current_state <- sample(1:num_states, 1, prob = P[current_state,])
  }
  
  # Ajustamos el número de registros a exactamente L
  chain <- chain[1:L]
  gl <- gl[1:L]
  time <- time[1:L]
  
  # Agregar ruido gaussiano a los valores de glucosa
  gl_noise <- pmax(gl + rnorm(L, mean = 0, sd = noise_sd), 0)
  
  # Formateamos el tiempo al formato solicitado
  time <- format(time, "%Y-%m-%d %H:%M:%S")
  
  # Crear un dataframe con las columnas solicitadas
  result <- data.frame(
    id = rep(id, L),
    time = time,
    gl = gl,
    gl_noise = gl_noise,  # Columna con el ruido agregado
    state = chain  # Usamos la secuencia generada por el modelo
  )
  
  return(result)
}
