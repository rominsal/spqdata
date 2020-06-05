#' @title A funcion to calculate Q CON PLOTS
#'
#' @description This function calculates Q, a measure of spatial association based on symbolic entropy.
#' @param fx vector de observaciones Factor o numérico?
#' @param mh m-historias
#' @param symb matrix con simbolos
#' @usage q_symb(fx, mh, symb)
#' @keywords spatial association, qualitative variable, symbolic entropy, symbols
#' @details Aquí Antonio escribe una linda historia
#' @return decir que cosas son las que devuelve
#'   \tabular{ll}{
#'     \code{pc} \tab Qp for combinations-totals symbols \cr
#'     \code{pc_pval} \tab  p-value for Qp \cr
#'     \code{qc} \tab Qc for combinations-totals symbols \cr
#'     \code{qc_pval} \tab  p-value for Qc \cr
#'     \code{p_symb} \tab  Matriz que lista los símbolos sin compactar \cr
#'     \code{efp_symb} \tab  Frecuancia absoluta de cada símbolo (no compactados)\cr
#'     \code{c_symb} \tab  Matriz que lista los símbolos COMPACTOS \cr
#'     \code{efc_symb} \tab  Frecuancia absoluta de cada símbolo (Compactado)\cr
#'     }
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Antonio  \tab \email{paez@@gmail.com} \cr
#'   Manolo  \tab \email{manuel.ruiz@@upct.es} \cr
#'   }
#'   @references
#'   \itemize{
#'     \item Ruiz, M., López, F., and Páez, A. (2010).
#'     Testing for spatial association of qualitative data using symbolic dynamics.
#'       \emph{Journal of Geographical Systems}, 12(3), 281-309.
#'     \item Ruiz, M., López, F., and Páez, A. (2010).
#'     Testing for spatial association of qualitative data using symbolic dynamics.
#'       \emph{Journal of Geographical Systems}, 12(3), 281-309.0.
#'   }
#' @seealso
#' \code{\link{dgp_spq}}, \code{\link{m_surr_no}}
#' @export
#' @examples
#'
#' # Load dataset
#' data("FastFood")
#' # Define coordinates
#' x <- cbind(FastFood.sf$Lat,FastFood.sf$Lon)
#' symb33 <- cr_symb(3,3) # Obtain symbols for 3 classes
#' # Obtain m-surroundings of size 3 (m=3), with degree of overlap of one (s=1)
#' m <- 3
#' s <- 1
#' mh31 <- m_surr_no(x,m,s)
#'
#' # Plot the m-sourronding
#' N <- dim(x)[1]
#' W <- matrix(0,ncol = N,nrow = N)
#' for (i in 1:dim(mh31)[1]){W[mh31[i,1],mh31[i,2:m]] <- 1}
#' g1 = igraph::graph.adjacency(W)
#' mylist <- list()
#' for (i in 1:dim(mh31)[1]){mylist[[i]] <- mh31[i,]}
#'
#' plot(g1,layout=layout.norm(x),edge.arrow.mode=0,edge.width=0,vertex.size=1,vertex.label="",vertex.color="red",vertex.label.font=0,mark.groups=mylist)
#'
#'
#'
#' # Obtain the Q test
#' W <- mat2listw(W)
#' results.31 <- q_symb_A(FastFood.sf$Type,mh31,symb33)

q_symb_A <- function(fx, mh, symb) {

  if (is.factor(fx)){
    levels(fx) <- as.character(1:length(levels(fx)))
    Y <- as.numeric(fx)
  }
  if (is.character(fx)){
    Y <- as.factor(fx)
    levels(fx) <- as.character(1:length(levels(fx)))
    Y <- as.numeric(fx)
  }
  if (is.numeric(fx)){
    Y <- fx
  }
  mdim = dim(mh)
  n1 <- mdim[1]  # n1: number of symbolized locations
  m <- mdim[2]  # m: size of m-surrounding
  N = length(Y)  # N: number of observations
  # Calculate the probabilities of each of K qualitative classes
  Ymax <- max(Y)
  pm <- rep(0, Ymax)  # Initialize probability vector
  for (k in 1:Ymax) {
    pm[k] <- sum(Y == k)/N  # Probabilities of each of K qualitative classes
  }

  np_symb <- nrow(symb$p_symb)  #np_symb: number of standard permutation symbols
  # Collapse the symbols to characters to then use as factors
  PSymb <- paste(symb$p_symb[, 1])
  for (i in 2:m) {
    PSymb <- paste(PSymb, symb$p_symb[, i])
  }
  # Standard symbols in equivalents by combination-totals
  p_symbi <- matrix(rep(0, np_symb * Ymax), ncol = Ymax)
  for (i in 1:np_symb) {
    for (j in 1:k) {
      p_symbi[i, j] = sum(symb$p_symb[i, ] == j)
    }
  }
  # Collapse the symbols to characters to then use as factors
  PSymbi <- paste(p_symbi[, 1])
  for (i in 2:Ymax) {
    PSymbi <- paste(PSymbi, p_symbi[, i])
  }

  nc_symb <- nrow(symb$c_symb)  # np_symb: number of reduced combination symbols
  # Collapse the symbols to characters to then use as factors
  CSymb <- paste(symb$c_symb[, 1])
  for (i in 2:Ymax) {
    CSymb <- paste(CSymb, symb$c_symb[, i])
  }
  n2 <- m
  Sp <- rep(0, n1)  # Initialize the vector of permutations symbols for all symbolized locations
  Sc <- rep(0, n1)  # Initialize the vector of combinations-totals symbols for all symbolized locations

  # Evaluate the probabilities of each symbol under the null, standard
  # permutation symbols.
  qp_symb <- rep(0, np_symb)
  for (k in 1:np_symb) {
    qp_symb[k] <- 1
    for (j in 1:m) {
      qp_symb[k] <- qp_symb[k] * pm[symb$p_symb[k, j]]  # frequency of symbols, is the joint probability of classes in symbols
    }
  }
  # Expected frequency of standard permutation symbols under the null
  sfp_0 <- n1 * qp_symb  # The probability of the symbol under the null times the number of symbolized locations
  # Symbolize locations
  Z <- matrix(Y[mh], ncol = m)
  for (i in 1:n1) {
    for (k in 1:np_symb) {
      if (all(symb$p_symb[k, ] == Z[i, ])) {
        Sp[i] = k  #S identifies, for each symbolized location, its corresponding symbol
      }
    }
  }

  # Empirical frequency of standard symbols
  efp_symb <- rep(0,np_symb)
  for (i in 1:np_symb){
    efp_symb[i] <- sum(Sp==i)
  }
  #efp_symb <- tabulate(Sp)

  ############################################################################
  # Plot 1
  ############################################################################
  # Standard error 95% for standard permutation symbols
  sep_symb <- qnorm(0.975, lower.tail = TRUE) * sqrt((sfp_0/n1 * (1 -
                                                                    sfp_0/n1))/n1)
  # Determine if empirical frequency is outside of intervals of
  # confidence
  sigp_symb <- (-1 * (efp_symb/n1 < (sfp_0/n1 - qnorm(0.975, lower.tail = TRUE) *
                                       sqrt((sfp_0/n1 * (1 - sfp_0/n1))/n1))) + 1 * (efp_symb/n1 > (sfp_0/n1 +
                                       qnorm(0.975, lower.tail = TRUE) * sqrt((sfp_0/n1 * (1 - sfp_0/n1))/n1))))
  sigp_symb <- factor(sigp_symb)
  if (any(sigp_symb==-1)){
    levels(sigp_symb)[levels(sigp_symb)=="-1"] <- "-"
  }
  if (any(sigp_symb==0)){
    levels(sigp_symb)[levels(sigp_symb)=="0"] <- "NS"
  }
  if (any(sigp_symb==1)){
    levels(sigp_symb)[levels(sigp_symb)=="1"] <- "+"
  }
  # Dataframe for plotting
  PSymb.df <- data.frame(PSymb, efp_symb/n1, sfp_0/n1, sep_symb, sigp_symb)
  # Create ggplot2 plot object
  p_symb_plot <- ggplot2::ggplot(PSymb.df) + ggplot2::geom_bar(ggplot2::aes(x = PSymb,
  y = efp_symb.n1, color = sigp_symb), stat = "identity") + ggplot2::geom_errorbar(ggplot2::aes(x = PSymb,
  y = sfp_0.n1, ymin = sfp_0.n1 - sep_symb, ymax = sfp_0.n1 + sep_symb,
  color = sigp_symb)) + ggplot2::labs(x = "Symbol (Permutations)",
  y = "Frequency", color = "Significance") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
  hjust = 0, size = 6), legend.position = "top") + ggplot2::theme_bw()

  # Evaluate the probabilities of each symbol under the null, combination
  # symbols.
  qc_symb <- rep(0, nc_symb)
  for (k in 1:nc_symb) {
    # qc_symb[k]=1 for(j in 1:Ymax){
    qc_symb[k] = sum((PSymbi == CSymb[k]) * qp_symb)  #frequency of combination symbols
    # }
  }
  # Expected frequency of standard permutation symbols under the null
  sfc_0 <- n1 * qc_symb  #The probability of the symbol under the null times the number of symbolized locations
  # Symbolize locations
  Zm <- matrix(rep(0, n1 * Ymax), ncol = Ymax)
  for (j in 1:n1) {
    for (k in 1:Ymax) {
      Zm[j, k] <- sum(Z[j, ] == k)
    }
  }
  for (i in 1:n1) {
    for (k in 1:nc_symb) {
      if (all(symb$c_symb[k, ] == Zm[i, ])) {
        Sc[i] = k  #S identifies, for each symbolized location, its corresponding symbol
      }
    }
  }

  # Empirical frequency of combination symbols
  efc_symb <- rep(0,nc_symb)
  for (i in 1:nc_symb){
    efc_symb[i] <- sum(Sc==i)
  }

  #efc_symb <- tabulate(Sc)

  ############################################################################
  # Standard error 95% for standard permutation symbols
  sec_symb <- qnorm(0.975, lower.tail = TRUE) * sqrt((sfc_0/n1 * (1 -
                                                                    sfc_0/n1))/n1)
  # Determine if empirical frequency is outside of intervals of
  # confidence
  sigc_symb <- (-1 * (efc_symb/n1 < (sfc_0/n1 - qnorm(0.975, lower.tail = TRUE) *
                                       sqrt((sfc_0/n1 * (1 - sfc_0/n1))/n1))) + 1 * (efc_symb/n1 > (sfc_0/n1 +
                                                                                                      qnorm(0.975, lower.tail = TRUE) * sqrt((sfc_0/n1 * (1 - sfc_0/n1))/n1))))
  sigc_symb <- factor(sigc_symb)
  if (any(sigc_symb==-1)){
    levels(sigc_symb)[levels(sigc_symb)=="-1"] <- "-"
  }
  if (any(sigc_symb==0)){
    levels(sigc_symb)[levels(sigc_symb)=="0"] <- "NS"
  }
  if (any(sigc_symb==1)){
    levels(sigc_symb)[levels(sigc_symb)=="1"] <- "+"
  }
  # Dataframe for plotting
  CSymb.df <- data.frame(CSymb, efc_symb/n1, sfc_0/n1, sec_symb, sigc_symb)
  # Create ggplot2 plot object
  c_symb_plot <- ggplot2::ggplot(CSymb.df) + ggplot2::geom_bar(ggplot2::aes(x = CSymb,
          y = efc_symb.n1, color = sigc_symb), stat = "identity") + ggplot2::geom_errorbar(ggplot2::aes(x = CSymb,
           y = sfc_0.n1, ymin = sfc_0.n1 - sec_symb, ymax = sfc_0.n1 + sec_symb,
         color = sigc_symb)) + ggplot2::labs(x = "Symbol (Combinations - Totals)",
          y = "Frequency", color = "Significance") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
           hjust = 0, size = 6), legend.position = "top") + ggplot2::theme_bw()

  # Calculate statistics With standard permutations symbols
  lnp <- rep(0, np_symb)
  for (i in 1:np_symb) {
    if (efp_symb[i] != 0) {
      lnp[i] = log(efp_symb[i]/n1)
    }
  }
  hmp = -sum((efp_symb/n1) * lnp)  #Empirical
  hmp_0 = sum((efp_symb/n1) * log(qp_symb))  #Expected under the null
  qp = -2 * n1 * sum(hmp + hmp_0)  #Likelihood ratio statistic for permutations symbols
  qp_pval = pchisq(qp, df = np_symb-1, lower.tail = FALSE)

  # With combinations-totals symbols
  lnc <- rep(0, nc_symb)
  for (i in 1:nc_symb) {
    if (efc_symb[i] != 0) {
      lnc[i] = log(efc_symb[i]/n1)
    }
  }
  hmc = -sum((efc_symb/n1) * lnc)  #Empirical
  hmc_0 = sum((efc_symb/n1) * log(qc_symb))  #Expected under the null
  qc = -2 * n1 * sum(hmc + hmc_0)  # Likelihood ratio statistic for combinations-totals symbols
  qc_pval = pchisq(qc, df = nc_symb-1, lower.tail = FALSE)

  # Return results
  results <- list(qp, qp_pval, qc, qc_pval, symb$p_symb,efp_symb,symb$c_symb,efc_symb, p_symb_plot, c_symb_plot)
  names(results)[1] <- "qp"
  names(results)[2] <- "qp_pval"
  names(results)[3] <- "qc"
  names(results)[4] <- "qc_pval"
  names(results)[5] <- "p_symb" # simbolos p (Sin compactar)
  names(results)[6] <- "efp_symb" # Distribución empírica de símbolos
  names(results)[7] <- "c_symb" # simbolos c (Compactos)
  names(results)[8] <- "efc_symb" # Distribución empírica de símbolos
  names(results)[9] <- "p_symb_plot" # simbolos c (Compactos)
  names(results)[10] <- "c_symb_plot" # Distribución empírica de símbolos

  return(results)  #To return a list of results a,b are two outputs

}
