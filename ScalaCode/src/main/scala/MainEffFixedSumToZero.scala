import java.io.File

/**
  * Saturated model for the main effects as fixed effects. Sum to 0 constraint added to solve the identifiability problem.
  */
object MainEffFixedSumToZero {

  def saturatedModelCatMainEffs(noOfIter: Int, y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], lalpha: Int, lbeta: Int, structure: DVStructure, ujPriorMean: Double, ujPriorVar: Double, nkPriorMean: Double, nkPriorVar:Double, muPrior: Double, varPrior: Double, a: Double, b:Double) = {

    val N = y.length

    val mat_mt = new DenseMatrix[Double](noOfIter, 2) //to store mu and tau
    val alphaCoefs = new DenseMatrix[Double](noOfIter, lalpha) //to store the coefficients for alpha
    val betaCoefs = new DenseMatrix[Double](noOfIter, lbeta) // to store the coefficients for beta

    val curAlpha = DenseVector.zeros[Double](lalpha) //the current values of the coefficients for every iteration. Initialised with 0.
    val curBeta = DenseVector.zeros[Double](lbeta) //the current values of the coefficients for every iteration. Initialised with 0.


    var mu = 0.0 //initialise the sampler at a=0
    var tau = 1.0 // initialise the sampler at tau=0

    //Initialise with 0 in the first row the matrices where we will store the values to start the Gibbs sampler from 0s
    mat_mt(0,::):= DenseVector(mu,tau).t
    alphaCoefs(0,::):= curAlpha.t
    betaCoefs(0,::):= curBeta.t
    val SYi= y.toArray.sum

    for(i<-1 until noOfIter ){

      //Update mu and tau
      //      mu = 3.10
      //      tau = 1.0
      val varMu= 1.0/((1.0/varPrior)+N*tau)
      val meanMu= ((muPrior/varPrior)+tau*(SYi - sumAllEffects(y,alpha, beta, curAlpha, curBeta, lalpha, lbeta)))/(1.0/varMu)
      mu= breeze.stats.distributions.Gaussian(meanMu,sqrt(varMu)).draw()
      tau= breeze.stats.distributions.Gamma(a+N/2.0, 1.0/(b+YminusEffects(y,alpha,beta,mu,curAlpha,curBeta)/2.0)).draw() //  !!!!TO SAMPLE FROM THE GAMMA DISTRIBUTION IN BREEZE THE β IS 1/β


      //        // Update alphaj
      for(jj <-0 until lalpha-1){
        val SYiuj= structure.calcAlphaSum(jj)
        val Nj= structure.calcAlphaLength(jj)
        val Sn= sumBetaEff(y,alpha,beta, jj, curBeta)
        val varPU = 1.0/((1.0/ujPriorVar)+ Nj*tau)
        val meanPU = ((ujPriorMean/ujPriorVar)+ tau*(SYiuj-Nj*mu-Sn))/(1.0/varPU)
        curAlpha(jj) = breeze.stats.distributions.Gaussian(meanPU, sqrt(varPU)).draw()
      }
      curAlpha(lalpha-1)= -sum(curAlpha(0 until lalpha-1))

      for(jk <-0 until lbeta-1){
        val SYink= structure.calcBetaSum(jk)
        val Nk= structure.calcBetaLength(jk)
        val Su= sumAlphaEff(y,alpha,beta, jk, curAlpha)
        val varPN = (1.0/((1.0/nkPriorVar)+ tau*Nk))
        val meanPN = ((nkPriorMean/nkPriorVar)+ tau*(SYink-Nk*mu-Su))/(1.0/varPN)
        curBeta(jk)= breeze.stats.distributions.Gaussian(meanPN, sqrt(varPN)).draw()
      }

      curBeta(lbeta-1)= -sum(curBeta(0 until lbeta-1))


      //curAlpha:= DenseVector(-1.0,4.0,-3.0)
      //curBeta:= DenseVector(-1.0,3.5,-2.0,-0.5)
      mat_mt(i,::):= DenseVector(mu,tau).t
      alphaCoefs(i,::):= curAlpha.t
      betaCoefs(i,::):= curBeta.t
    }
    (mat_mt,alphaCoefs,betaCoefs)
  }
  /**
    * Add all the beta effects for a given alpha.
    */
  def sumBetaEff(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], alphaIndex: Int, betaEff: DenseVector[Double]):Double={
    var sum=0.0
    for(i <- 0 until y.length){
      if(alpha(i)== alphaIndex+1){
        sum= sum + betaEff(beta(i)-1)
      }
    }
    sum
  }
  /**
    * Add all the alpha effects for a given beta.
    */
  def sumAlphaEff(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], betaIndex: Int, alphaEff: DenseVector[Double]):Double={
    var sum=0.0
    for(i <- 0 until y.length){
      if(beta(i)== betaIndex+1){
        sum= sum + alphaEff(alpha(i)-1)
      }
    }
    sum
  }

  /**
    * Calculate the Yi-mu-u_eff-n_eff. To be used in estimating tau
    */
  def YminusEffects(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], mu: Double, alphaEff: DenseVector[Double],betaEff: DenseVector[Double]):Double={
    val N= y.length
    var YminusMuBetaAlphaEff = DenseVector.zeros[Double](N)
    for (i <- 0 until N){
      YminusMuBetaAlphaEff(i)= y(i)-alphaEff(alpha(i)-1)-betaEff(beta(i)-1)-mu
    }
    sqr(YminusMuBetaAlphaEff).toArray.sum
  }

  /**
    * Calculate the sum of all the alpha and all the beta effects for all the observations.
    */
  def sumAllEffects(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int], alphaEff: DenseVector[Double], betaEff: DenseVector[Double], lalpha: Int, lbeta: Int):Double={
    var sumBeta=0.0
    var sumAlpha=0.0
    // For alpha effects
    for(i <- 0 until lbeta){
      sumAlpha= sumAlpha + sumAlphaEff(y,alpha,beta, i, alphaEff)
    }
    sum
    // For beta effects
    for(i <- 0 until lalpha){
      sumBeta= sumBeta + sumBetaEff(y,alpha,beta, i, betaEff)
    }
    sumAlpha+sumBeta
  }

  def sqr(x: DenseVector[Double]) = x * x

  /**
    * SATURATED MODEL. MAIN EFFECTS AS FIXED. SUM TO 0.
    * Saturated model main effects as fixed with sum to Zero constraint
    */
  def main(args: Array[String]): Unit = {
    // Read the data
    val data = csvread(new File("./simulNoInter.csv"))

    val sampleSize = data.rows
    val y = data(::, 0)
    val alpha = data(::, 1).map(_.toInt)
    val beta = data(::, 2).map(_.toInt)
    val structure = new DVStructure(y, alpha, beta)

    // Parameters
    val lalpha = data(::, 1).toArray.distinct.length
    val lbeta = data(::, 2).toArray.distinct.length
    val noOfIters = 100000
    val aPrior = 1
    val bPrior = 0.0001
    val ujPriorMean = 0.0
    val ujPriorVar = 1 / 0.0001
    val nkPriorMean = 0.0
    val nkPriorVar = 1 / 0.0001
    val muPrior = 0.0
    val varPrior = 1.0 / 0.0001
    val a = 1
    val b = 0.0001


    val (test_mtS0,alpha_estS0,beta_estS0)= saturatedModelCatMainEffs(noOfIters, y, alpha, beta, lalpha, lbeta, structure, ujPriorMean, ujPriorVar, nkPriorMean, nkPriorVar, muPrior, varPrior, a, b)
    val mtS0= mean(test_mtS0(::,*)).t
    val alphaEstimS0= mean(alpha_estS0(::,*)).t
    val betaEstimS0= mean(beta_estS0(::,*)).t

    println("Results: ")
    println("mu and tau:" + mtS0)
    println("Estimates for alpha: "+ alphaEstimS0)
    println("Estimates for beta: "+ betaEstimS0)

    //Plot Results
    val resMatS0= DenseMatrix.horzcat(test_mtS0,alpha_estS0, beta_estS0)
    val plotS0= new PlotResults
    plotS0.plotResults(resMatS0, List("Mean","tau", "alpha 1", "alpha 2", "alpha 3","beta 1","beta 2","beta 3","beta 4"))
  }


}
