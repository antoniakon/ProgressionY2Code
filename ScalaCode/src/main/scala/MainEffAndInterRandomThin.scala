import java.io.File

/**
  * Created by Antonia Kontaratou.
  * Saturated model for Main + interaction effects as random effects by putting a gamma prior on the precision of the effects. Thinning factor.
  * This version uses the new notation with alpha and beta to be consistent with the maths.
  */
object MainEffAndInterRandomThin {

  //---------------------- Main effects ---------------------------//
  // the model is: Xijk|mu,aj,bk,gjk, tau~N(mu+aj+bk+gjk,tau^-1)

  def satModelCatMainEfAndInterThin(noOfIter: Int, thin:Int, y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], nj: Int, nk: Int, structure: DVStructure, alphaPriorMean: Double, alphaPriorVar: Double, betaPriorMean: Double, betaPriorVar: Double, mu0: Double, tau0: Double, a: Double, b:Double, gammaPriorMean: Double, aPrior: Double, bPrior:Double) = {

    val N = y.length //the number of observations
    val sampleNo= noOfIter/thin+1 //the number of samples created from the MCMC
    val njk= nj*nk // the number of levels of interactions

    val mat_mt = new DenseMatrix[Double](sampleNo, 2) //to store mu and tau
    val alphaCoefs = new DenseMatrix[Double](sampleNo, nj) //to store the coefficients for variable alpha
    val betaCoefs = new DenseMatrix[Double](sampleNo, nk) // to store the coefficients for variable beta
    val gammaCoefs = new DenseMatrix[Double](sampleNo, njk) // to store the coefficients for the interactions, variable gamma
    val SumX = y.toArray.sum // the sum of the values of all the observations

    val curAlpha = DenseVector.zeros[Double](nj) //the current values of the coefficients for every iteration. Initialised with 0.
    val curBeta = DenseVector.zeros[Double](nk)
    val curGamma = DenseMatrix.zeros[Double](nj, nk)

    var mu = 0.0 //initialise the sampler at mu=0
    var tau = 1.0 // initialise the sampler at tau=1.0

    var tauAlpha= 1.0 //initialise the precision for alpha (this is actually tau_alpha)
    var tauBeta= 1.0 //initialise the precision for beta (this is actually tau_beta)
    var tauGamma= 1.0 //initialise the precision for beta (this is actually tau_gamma)


    var sumaj=0.0 //this will be used for sampling the precision of alpha from the FCD. Here it is initialised at 0.
    var sumbk=0.0 //this will be used for sampling the precision of beta from the FCD. Here it is initialised at 0.
    var sumgjk=0.0 //this will be used for sampling the precision of the interactions gamma from the FCD. Here it is initialised at 0.

    var ind= 0 //index used for thinning

    //Initialise with 0 in the first row the matrices where we will store the values to start the Gibbs sampler from 0s
    mat_mt(0, ::) := DenseVector(mu, tau).t
    alphaCoefs(0, ::) := curAlpha.t
    betaCoefs(0, ::) := curBeta.t
    gammaCoefs(0, ::) := curGamma.t.toDenseVector.t

    for (i <- 1 until noOfIter) {
      println("iter :" + i)

      for (j<-0 until nj){
        sumaj= sumaj+pow((curAlpha(j)-alphaPriorMean),2) //estimate the sum used in sampling from Gamma distribution for the precision of alpha
      }

      for (k<-0 until nk){
        sumbk= sumbk+pow((curBeta(k)-betaPriorMean),2) //estimate the sum used in sampling from Gamma distribution for the precision of beta
      }

      for (j<-0 until nj) {
        for (k <- 0 until nk) {
          sumgjk = sumgjk + pow((curGamma(j,k) - gammaPriorMean), 2) //estimate the sum used in sampling from Gamma distribution for the precision of gamma/interacions
        }
      }

      tauAlpha= breeze.stats.distributions.Gamma(aPrior+nj/2.0, 1.0/(bPrior+0.5*sumaj)).draw() //sample the precision of alpha from gamma
      tauBeta= breeze.stats.distributions.Gamma(aPrior+nk/2.0, 1.0/(bPrior+0.5*sumbk)).draw() // sample the precision of beta from gamma
      tauGamma= breeze.stats.distributions.Gamma(aPrior+njk/2.0, 1.0/(bPrior+0.5*sumgjk)).draw() // sample the precision of he interactions gamma from gamma Distribition

      //Update mu and tau
      //mu = 3.10
      //tau = 1.0
      val varMu = 1.0/(tau0+N*tau) //the variance for mu
      val meanMu = (mu0*tau0 + tau*(SumX - sumAllMainInterEff(y, alpha, beta, curAlpha, curBeta, nj, nk, curGamma)))*varMu
      mu = breeze.stats.distributions.Gaussian(meanMu, sqrt(varMu)).draw()
      tau = breeze.stats.distributions.Gamma(a+N/2.0, 1.0 / (b +0.5*YminusMuAndEffects(y, alpha, beta, mu, curAlpha, curBeta, curGamma))).draw() //  !!!!TO SAMPLE FROM THE GAMMA DISTRIBUTION IN BREEZE THE β IS 1/β

      // Update alphaj
      for (j <- 0 until nj) {
        val SXalphaj = structure.calcAlphaSum(j) // the sum of the observations that have alpha==j
        val Nj = structure.calcAlphaLength(j) // the number of the observations that have alpha==j
        val SumBeta = sumBetaEffGivenAlpha(y, alpha, beta, j, curBeta) //the sum of the beta effects given alpha
        val SinterAlpha = sumInterEffGivenAlpha(y, alpha, beta, j, curGamma) //the sum of the gamma/interaction effects given alpha
        val varPalpha = 1.0 / (tauAlpha + tau * Nj) //the variance for alphaj
        val meanPalpha = (alphaPriorMean*tauAlpha + tau * (SXalphaj - Nj * mu - SumBeta - SinterAlpha)) *varPalpha //the mean for alphaj
        curAlpha(j) = breeze.stats.distributions.Gaussian(meanPalpha, sqrt(varPalpha)).draw()
      }

      // Update betak
      for (k <- 0 until nk) {
        val SXbetak = structure.calcBetaSum(k) // the sum of the observations that have beta==k
        val Nk = structure.calcBetaLength(k) // the number of the observations that have beta==k
        val SumAlpha = sumAlphaGivenBeta(y, alpha, beta, k, curAlpha) //the sum of the alpha effects given beta
        val SinterBeta = sumInterEffGivenBeta(y, alpha, beta, k, curGamma) //the sum of the gamma/interaction effects given beta
        val varPbeta = 1.0 / (tauBeta + tau * Nk) //the variance for betak
        val meanPbeta = (betaPriorMean*tauBeta + tau * (SXbetak - Nk * mu - SumAlpha - SinterBeta))*varPbeta //the mean for betak
        curBeta(k) = breeze.stats.distributions.Gaussian(meanPbeta, sqrt(varPbeta)).draw()
      }

      // Update Interaction terms
      for (j <- 0 until nj) {
        for (k <- 0 until nk) {
          val Njk= structure.getDVList(j,k).length // the number of the observations that have alpha==j and beta==k
          val SXjk = structure.getDVList(j,k).sum // the sum of the observations that have alpha==j and beta==k
          val varPInter = 1.0 / (tauGamma + tau * Njk) //the variance for gammajk
          val meanPInter=(gammaPriorMean*tauGamma + tau * (SXjk - Njk*(mu+curAlpha(j)+curBeta(k))))*varPInter
          curGamma(j,k) = breeze.stats.distributions.Gaussian(meanPInter, sqrt(varPInter)).draw()
        }
      }

      if((i%thin).equals(0)){
        ind=ind+1
        mat_mt(ind, ::) := DenseVector(mu, tau).t
        alphaCoefs(ind, ::) := curAlpha.t
        betaCoefs(ind, ::) := curBeta.t
        gammaCoefs(ind, ::) := curGamma.t.toDenseVector.t
      }
//           curAlpha:= DenseVector(-1.0,4.0,-3.0)
//           curBeta:= DenseVector(-1.0,3.5,-2.0,-0.5)
      //curGamma:= DenseMatrix((-0.5,-2.0,0.0,2.5),  (0.0,-0.5,-1.0,1.5), (0.5,2.5,1.0,-4.0))
      sumaj=0.0
      sumbk=0.0
      sumgjk=0.0
    }
    (mat_mt, alphaCoefs, betaCoefs, gammaCoefs)
  }

  /**
    * Add all the beta effects for a given alpha.
    */
  def sumBetaEffGivenAlpha(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], alphaIndex: Int, betaEff: DenseVector[Double]):Double={
    val N= y.length
    var sum=0.0
    for(i <- 0 until N){                  // run through all the observations
      if(alpha(i)== alphaIndex+1){        // if alpha of the current observation == current alphaIndex [+1 because of the difference in the dataset notation (alpha=1,2,...) and Scala indexing that starts from 0]
        sum= sum + betaEff(beta(i)-1)     // add to the sum the current effect of the observation's beta (-1 for consistency with the Scala vector indexing again)
      }
    }
    sum
  }
  /**
    * Add all the alpha effects for a given beta.
    */
  def sumAlphaGivenBeta(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], betaIndex: Int, alphaEff: DenseVector[Double]):Double={
    val N= y.length
    var sum=0.0
    for(i <- 0 until N){
      if(beta(i)== betaIndex+1){
        sum= sum + alphaEff((alpha(i)-1))
      }
    }
    sum
  }
  /**
    * Calculate the sum of all the alpha and all the beta effects for all the observations.
    */
  def sumAllMainInterEff(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], alphaEff: DenseVector[Double], betaEff: DenseVector[Double], nj: Int, nk: Int, interEff: DenseMatrix[Double]):Double={
    var sumBeta=0.0
    var sumAlpha=0.0
    var sumInter=0.0
    // For alpha effects
    for(i <- 0 until nk){ //through all beta effects
      sumAlpha= sumAlpha + sumAlphaGivenBeta(y,alpha,beta, i, alphaEff)
    }
    sum
    // For beta effects
    for(i <- 0 until nj){ //through all alpha effects
      sumBeta= sumBeta + sumBetaEffGivenAlpha(y,alpha,beta, i, betaEff)
    }
    // For Interaction effects
    for(i <- 0 until nj){
      for(j <- 0 until nk){
        sumInter= sumInter + sumInterEff(y,alpha,beta, i, j, interEff)
      }
    }

    sumAlpha+sumBeta+sumInter
  }

  /**
    * Add all the interaction effects for a given alpha and a given beta.
    */
  def sumInterEff(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], alphaIndex: Int, betaIndex: Int, interEff: DenseMatrix[Double]):Double={
    var sum=0.0
    for(i <- 0 until y.length){
      if((alpha(i)== alphaIndex+1)&(beta(i)== betaIndex+1)){
        sum = sum + interEff((alpha(i)-1), (beta(i)-1))
      }
    }
    sum
  }

  /**
    * Add all the interaction effects for a given alpha.
    */
  def sumInterEffGivenAlpha(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], alphaIndex: Int, interEff: DenseMatrix[Double]):Double={
    val N= y.length
    var sum=0.0
    for(i <- 0 until N){
      if((alpha(i)== alphaIndex+1)){
        sum = sum + interEff((alpha(i)-1), (beta(i)-1))
      }
    }
    sum
  }
  /**
    * Add all the interaction effects for a given beta.
    */
  def sumInterEffGivenBeta(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], betaIndex: Int, interEff: DenseMatrix[Double]):Double={
    var sum=0.0
    for(i <- 0 until y.length){
      if((beta(i)== betaIndex+1)){
        sum = sum + interEff((alpha(i)-1), (beta(i)-1))
      }
    }
    sum
  }
  /**
    * Calculate the Yi-mu-u_eff-n_eff- inter_effe. To be used in estimating tau
    */
  def YminusMuAndEffects(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], mu: Double, alphaEff: DenseVector[Double], betaEff: DenseVector[Double], interEff: DenseMatrix[Double]):Double={
    val N= y.length
    val YminusMuNatUniEff = DenseVector.zeros[Double](N)
    for (i <- 0 until N){
      YminusMuNatUniEff(i)= y(i)-mu-alphaEff(alpha(i)-1)-betaEff(beta(i)-1)-interEff(alpha(i)-1, beta(i)-1)
    }
    sqr(YminusMuNatUniEff).toArray.sum
  }

  def sqr(x: DenseVector[Double]) = x * x

  // Calculation of the execution time
  def time[A](f: => A) = {
    val s = System.nanoTime
    val ret = f
    println( "time: " +(System.nanoTime-s)/1e6+ "ms")
    ret
  }

  def main(args: Array[String]): Unit = {
    // Read the data
    val data = csvread(new File("./simulInter.csv"))
    val sampleSize = data.rows
    val y = data(::, 0)
    val alpha = data(::, 1).map(_.toInt)
    val beta = data(::, 2).map(_.toInt)
    val structure = new DVStructure(y, alpha, beta)

    // Parameters
    val nj = data(::, 1).toArray.distinct.length
    val nk = data(::, 2).toArray.distinct.length
    val noOfIters = 1000000
    val thin = 100
    val aPrior = 1
    val bPrior = 0.0001
    val alphaPriorMean = 0.0
    val alphaPriorTau = 0.0001
    val betaPriorMean = 0.0
    val betaPriorTau = 0.0001
    val mu0 = 0.0
    val tau0 = 0.0001
    val a = 1
    val b = 0.0001
    val interPriorMean = 0.0 //common mean for all the interaction effects

    val (test_mtInter, alpha_estInter, beta_estInter, inter_est) = time(satModelCatMainEfAndInterThin(noOfIters, thin, y, alpha, beta, nj, nk, structure, alphaPriorMean, alphaPriorTau, betaPriorMean, betaPriorTau, mu0, tau0, a, b, interPriorMean, aPrior, bPrior))
    val mt = mean(test_mtInter(::, *)).t
    val alphaEstim = mean(alpha_estInter(::, *)).t
    val betaEstim = mean(beta_estInter(::, *)).t
    val interEstim = mean(inter_est(::, *)).t

    // Save the results to a csv file
    val mergedMatrix = DenseMatrix.horzcat(test_mtInter, alpha_estInter, beta_estInter, inter_est)
    val outputFIle = new File("./mainAndInterRandomThin.csv")
    breeze.linalg.csvwrite(outputFIle, mergedMatrix, separator = ',')

    println("Results: ")
    println("mu and tau:" + mt)
    println("Estimates for alpha: " + alphaEstim)
    println("Estimates for beta: " + betaEstim)
    println("Estimates for interactions: " + interEstim)

    //Plot Results
    val resMatR = DenseMatrix.horzcat(test_mtInter, alpha_estInter, beta_estInter, inter_est)
    val plotR = new PlotResults
    plotR.plotResults(resMatR, List("Mean", "tau", "alpha 1", "alpha 2", "alpha 3", "beta 1", "beta 2", "beta 3", "beta 4", "a1-b1", "a1-b2", "a1-b3", "a4-b4", "a2-b1", "a2-b2", "a2-b3", "a2-b4", "a3-b1", "a3-b2", "a3-b3", "a3-b4"))
  }

}
